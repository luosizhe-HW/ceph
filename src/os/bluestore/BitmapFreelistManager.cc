// -*- mode:C++; tab-width:8; c-basic-offset:2; indent-tabs-mode:t -*-
// vim: ts=8 sw=2 smarttab

#include "BitmapFreelistManager.h"
#include "kv/KeyValueDB.h"
#include "os/kv.h"

#include "common/debug.h"

#define dout_context cct
#define dout_subsys ceph_subsys_bluestore
#undef dout_prefix
#define dout_prefix *_dout << "freelist "

void make_offset_key(uint64_t offset, std::string *key)
{
  key->reserve(10);
  _key_encode_u64(offset, key);
}

struct XorMergeOperator : public KeyValueDB::MergeOperator {
  void merge_nonexistent(
    const char *rdata, size_t rlen, std::string *new_value) override {
    *new_value = std::string(rdata, rlen);
  }
  void merge(
    const char *ldata, size_t llen,
    const char *rdata, size_t rlen,
    std::string *new_value) override {
    ceph_assert(llen == rlen);
    *new_value = std::string(ldata, llen);
    for (size_t i = 0; i < rlen; ++i) {
      (*new_value)[i] ^= rdata[i];
    }
  }
  // We use each operator name and each prefix to construct the
  // overall RocksDB operator name for consistency check at open time.
  const char *name() const override {
    return "bitwise_xor";
  }
};

void BitmapFreelistManager::setup_merge_operator(KeyValueDB *db, string prefix)
{
  std::shared_ptr<XorMergeOperator> merge_op(new XorMergeOperator);
  db->set_merge_operator(prefix, merge_op);
}

BitmapFreelistManager::BitmapFreelistManager(CephContext* cct,
					     string meta_prefix,
					     string bitmap_prefix)
  : FreelistManager(cct),
    meta_prefix(meta_prefix),
    bitmap_prefix(bitmap_prefix),
    enumerate_bl_pos(0)
{
}

int BitmapFreelistManager::create(uint64_t new_size, uint64_t granularity,
				  KeyValueDB::Transaction txn)
{
  bytes_per_block = granularity;
  ceph_assert(isp2(bytes_per_block));
  size = p2align(new_size, bytes_per_block);
  blocks_per_key = cct->_conf->bluestore_freelist_blocks_per_key;

  _init_misc();

  blocks = size / bytes_per_block;
  if (blocks / blocks_per_key * blocks_per_key != blocks) {
    blocks = (blocks / blocks_per_key + 1) * blocks_per_key;
    dout(10) << __func__ << " rounding blocks up from 0x" << std::hex << size
	     << " to 0x" << (blocks * bytes_per_block)
	     << " (0x" << blocks << " blocks)" << std::dec << dendl;
    // set past-eof blocks as allocated
    _xor(size, blocks * bytes_per_block - size, txn);
  }
  dout(10) << __func__
	   << " size 0x" << std::hex << size
	   << " bytes_per_block 0x" << bytes_per_block
	   << " blocks 0x" << blocks
	   << " blocks_per_key 0x" << blocks_per_key
	   << std::dec << dendl;
  {
    bufferlist bl;
    encode(bytes_per_block, bl);
    txn->set(meta_prefix, "bytes_per_block", bl);
  }
  {
    bufferlist bl;
    encode(blocks_per_key, bl);
    txn->set(meta_prefix, "blocks_per_key", bl);
  }
  {
    bufferlist bl;
    encode(blocks, bl);
    txn->set(meta_prefix, "blocks", bl);
  }
  {
    bufferlist bl;
    encode(size, bl);
    txn->set(meta_prefix, "size", bl);
  }
  return 0;
}

int BitmapFreelistManager::expand(uint64_t new_size, KeyValueDB::Transaction txn)
{
  assert(new_size > size);
  ceph_assert(isp2(bytes_per_block));

  uint64_t blocks0 = size / bytes_per_block;
  if (blocks0 / blocks_per_key * blocks_per_key != blocks0) {
    blocks0 = (blocks0 / blocks_per_key + 1) * blocks_per_key;
    dout(10) << __func__ << " rounding blocks up from 0x" << std::hex << size
	     << " to 0x" << (blocks0 * bytes_per_block)
	     << " (0x" << blocks0 << " blocks)" << std::dec << dendl;
    // reset past-eof blocks to unallocated
    _xor(size, blocks0 * bytes_per_block - size, txn);
  }

  size = p2align(new_size, bytes_per_block);
  blocks = size / bytes_per_block;

  if (blocks / blocks_per_key * blocks_per_key != blocks) {
    blocks = (blocks / blocks_per_key + 1) * blocks_per_key;
    dout(10) << __func__ << " rounding blocks up from 0x" << std::hex << size
	     << " to 0x" << (blocks * bytes_per_block)
	     << " (0x" << blocks << " blocks)" << std::dec << dendl;
    // set past-eof blocks as allocated
    _xor(size, blocks * bytes_per_block - size, txn);
  }

  dout(10) << __func__
	   << " size 0x" << std::hex << size
	   << " bytes_per_block 0x" << bytes_per_block
	   << " blocks 0x" << blocks
	   << " blocks_per_key 0x" << blocks_per_key
	   << std::dec << dendl;
  {
    bufferlist bl;
    encode(blocks, bl);
    txn->set(meta_prefix, "blocks", bl);
  }
  {
    bufferlist bl;
    encode(size, bl);
    txn->set(meta_prefix, "size", bl);
  }
  return 0;
}

int BitmapFreelistManager::init(KeyValueDB *kvdb)
{
  dout(1) << __func__ << dendl;

  KeyValueDB::Iterator it = kvdb->get_iterator(meta_prefix);
  it->lower_bound(string());

  // load meta
  while (it->valid()) {
    string k = it->key();
    if (k == "bytes_per_block") {
      bufferlist bl = it->value();
      auto p = bl.cbegin();
      decode(bytes_per_block, p);
      dout(10) << __func__ << " bytes_per_block 0x" << std::hex
	       << bytes_per_block << std::dec << dendl;
    } else if (k == "blocks") {
      bufferlist bl = it->value();
      auto p = bl.cbegin();
      decode(blocks, p);
      dout(10) << __func__ << " blocks 0x" << std::hex << blocks << std::dec
	       << dendl;
    } else if (k == "size") {
      bufferlist bl = it->value();
      auto p = bl.cbegin();
      decode(size, p);
      dout(10) << __func__ << " size 0x" << std::hex << size << std::dec
	       << dendl;
    } else if (k == "blocks_per_key") {
      bufferlist bl = it->value();
      auto p = bl.cbegin();
      decode(blocks_per_key, p);
      dout(10) << __func__ << " blocks_per_key 0x" << std::hex << blocks_per_key
	       << std::dec << dendl;
    } else {
      derr << __func__ << " unrecognized meta " << k << dendl;
      return -EIO;
    }
    it->next();
  }

  dout(10) << __func__ << std::hex
	   << " size 0x" << size
	   << " bytes_per_block 0x" << bytes_per_block
	   << " blocks 0x" << blocks
	   << " blocks_per_key 0x" << blocks_per_key
	   << std::dec << dendl;
  _init_misc();
  return 0;
}

void BitmapFreelistManager::_init_misc()
{
  bufferptr z(blocks_per_key >> 3);
  memset(z.c_str(), 0xff, z.length());
  all_set_bl.clear();
  all_set_bl.append(z);

  block_mask = ~(bytes_per_block - 1);

  bytes_per_key = bytes_per_block * blocks_per_key;
  key_mask = ~(bytes_per_key - 1);
  dout(10) << __func__ << std::hex << " bytes_per_key 0x" << bytes_per_key
	   << ", key_mask 0x" << key_mask << std::dec
	   << dendl;
}

void BitmapFreelistManager::shutdown()
{
  dout(1) << __func__ << dendl;
}

void BitmapFreelistManager::enumerate_reset()
{
  std::lock_guard l(lock);
  enumerate_offset = 0;
  enumerate_bl_pos = 0;
  enumerate_bl.clear();
  enumerate_p.reset();
}

int get_next_clear_bit(bufferlist& bl, int start)
{
  const char *p = bl.c_str();
  int bits = bl.length() << 3;
  while (start < bits) {
    // byte = start / 8 (or start >> 3)
    // bit = start % 8 (or start & 7)
    unsigned char byte_mask = 1 << (start & 7);
    if ((p[start >> 3] & byte_mask) == 0) {
      return start;
    }
    ++start;
  }
  return -1; // not found
}

int get_next_set_bit(bufferlist& bl, int start)
{
  const char *p = bl.c_str();
  int bits = bl.length() << 3;
  while (start < bits) {
    int which_byte = start / 8;
    int which_bit = start % 8;
    unsigned char byte_mask = 1 << which_bit;
    if (p[which_byte] & byte_mask) {
      return start;
    }
    ++start;
  }
  return -1; // not found
}

bool BitmapFreelistManager::enumerate_next(KeyValueDB *kvdb, uint64_t *offset, uint64_t *length)
{
  std::lock_guard l(lock);

  // initial base case is a bit awkward
  if (enumerate_offset == 0 && enumerate_bl_pos == 0) {
    dout(10) << __func__ << " start" << dendl;
    enumerate_p = kvdb->get_iterator(bitmap_prefix);
    enumerate_p->lower_bound(string());
    // we assert that the first block is always allocated; it's true,
    // and it simplifies our lives a bit.
    ceph_assert(enumerate_p->valid());
    string k = enumerate_p->key();
    const char *p = k.c_str();
    _key_decode_u64(p, &enumerate_offset);
    enumerate_bl = enumerate_p->value();
    ceph_assert(enumerate_offset == 0);
    ceph_assert(get_next_set_bit(enumerate_bl, 0) == 0);
  }

  if (enumerate_offset >= size) {
    dout(10) << __func__ << " end" << dendl;
    return false;
  }

  // skip set bits to find offset
  while (true) {
    enumerate_bl_pos = get_next_clear_bit(enumerate_bl, enumerate_bl_pos);
    if (enumerate_bl_pos >= 0) {
      *offset = _get_offset(enumerate_offset, enumerate_bl_pos);
      dout(30) << __func__ << " found clear bit, key 0x" << std::hex
	       << enumerate_offset << " bit 0x" << enumerate_bl_pos
	       << " offset 0x" << *offset
	       << std::dec << dendl;
      break;
    }
    dout(30) << " no more clear bits in 0x" << std::hex << enumerate_offset
	     << std::dec << dendl;
    enumerate_p->next();
    enumerate_bl.clear();
    if (!enumerate_p->valid()) {
      enumerate_offset += bytes_per_key;
      enumerate_bl_pos = 0;
      *offset = _get_offset(enumerate_offset, enumerate_bl_pos);
      break;
    }
    string k = enumerate_p->key();
    dout(30) << "==k==:" << k << "k length:" << k.length() << dendl;
    //================================================================
    int flag = 0; 
    while(k.length()==10)
    {
      dout(30) << "01." << dendl;
      enumerate_p->next();
      enumerate_bl.clear();
      if (!enumerate_p->valid()){
        dout(30) << "02." << dendl;
        enumerate_offset += bytes_per_key;
        enumerate_bl_pos = 0;
        *offset = _get_offset(enumerate_offset, enumerate_bl_pos);
        flag = 1;
        break;
      }
      else{
        dout(30) << "03." << dendl;
        k = enumerate_p->key();
        if(k.length()==10)
          continue;
        else{
          dout(30) << "04." << dendl;
          flag = 2;
          break;
        }
      }
    }
    if(flag==1)
      break;
    //===========================================================
    const char *p = k.c_str();
    uint64_t next = enumerate_offset + bytes_per_key;
    _key_decode_u64(p, &enumerate_offset);
    enumerate_bl = enumerate_p->value();
    enumerate_bl_pos = 0;
    if (enumerate_offset > next) {
      dout(30) << " no key at 0x" << std::hex << next << ", got 0x"
	       << enumerate_offset << std::dec << dendl;
      *offset = next;
      break;
    }
  }

  // skip clear bits to find the end
  uint64_t end = 0;
  if (enumerate_p->valid()) {
    while (true) {
      enumerate_bl_pos = get_next_set_bit(enumerate_bl, enumerate_bl_pos);
      if (enumerate_bl_pos >= 0) {
	end = _get_offset(enumerate_offset, enumerate_bl_pos);
	dout(30) << __func__ << " found set bit, key 0x" << std::hex
		 << enumerate_offset << " bit 0x" << enumerate_bl_pos
		 << " offset 0x" << end << std::dec << dendl;
	end = std::min(get_alloc_units() * bytes_per_block, end);
	*length = end - *offset;
        dout(10) << __func__ << std::hex << " 0x" << *offset << "~" << *length
		 << std::dec << dendl;
	return true;
      }
      dout(30) << " no more set bits in 0x" << std::hex << enumerate_offset
	       << std::dec << dendl;
      enumerate_p->next();
      enumerate_bl.clear();
      enumerate_bl_pos = 0;
      if (!enumerate_p->valid()) {
	break;
      }
      string k = enumerate_p->key();
      //===================================================================
      int flag = 0;
      while(k.length()==10)
      {
        dout(30) << "001." << dendl;
        enumerate_p->next();
        enumerate_bl.clear();
        enumerate_bl_pos = 0;
        if(!enumerate_p->valid()){
          dout(30) << "002." << dendl;
          flag = 1;
          break;
        }
        else{
          dout(30) << "003." << dendl;
          k = enumerate_p->key();
          if(k.length()==10)
            continue;
          else{
            dout(30) << "004." << dendl;
            flag=2;
            break;
          }
        }
      }
      if(flag==1)
        break;
      //========================================================

      const char *p = k.c_str();
      _key_decode_u64(p, &enumerate_offset);
      enumerate_bl = enumerate_p->value();
    }
  }

  if (enumerate_offset < size) {
    end = get_alloc_units() * bytes_per_block;
    *length = end - *offset;
    dout(10) << __func__ << std::hex << " 0x" << *offset << "~" << *length
	     << std::dec << dendl;
    enumerate_offset = size;
    enumerate_bl_pos = blocks_per_key;
    return true;
  }

  dout(10) << __func__ << " end" << dendl;
  return false;
}

void BitmapFreelistManager::dump(KeyValueDB *kvdb)
{
  enumerate_reset();
  uint64_t offset, length;
  while (enumerate_next(kvdb, &offset, &length)) {
    dout(20) << __func__ << " 0x" << std::hex << offset << "~" << length
	     << std::dec << dendl;
  }
}

void BitmapFreelistManager::_verify_range(KeyValueDB *kvdb,
					  uint64_t offset, uint64_t length,
					  int val)
{
  unsigned errors = 0;
  uint64_t first_key = offset & key_mask;
  uint64_t last_key = (offset + length - 1) & key_mask;
  if (first_key == last_key) {
    string k;
    make_offset_key(first_key, &k);
    bufferlist bl;
    kvdb->get(bitmap_prefix, k, &bl);
    if (bl.length() > 0) {
      const char *p = bl.c_str();
      unsigned s = (offset & ~key_mask) / bytes_per_block;
      unsigned e = ((offset + length - 1) & ~key_mask) / bytes_per_block;
      for (unsigned i = s; i <= e; ++i) {
	int has = !!(p[i >> 3] & (1ull << (i & 7)));
	if (has != val) {
	  derr << __func__ << " key 0x" << std::hex << first_key << " bit 0x"
	       << i << " has 0x" << has << " expected 0x" << val
	       << std::dec << dendl;
	  ++errors;
	}
      }
    } else {
      if (val) {
	derr << __func__ << " key 0x" << std::hex << first_key
	     << " not present, expected 0x" << val << std::dec << dendl;
	++errors;
      }
    }
  } else {
    // first key
    {
      string k;
      make_offset_key(first_key, &k);
      bufferlist bl;
      kvdb->get(bitmap_prefix, k, &bl);
      if (bl.length()) {
	const char *p = bl.c_str();
	unsigned s = (offset & ~key_mask) / bytes_per_block;
	unsigned e = blocks_per_key;
	for (unsigned i = s; i < e; ++i) {
	  int has = !!(p[i >> 3] & (1ull << (i & 7)));
	  if (has != val) {
	    derr << __func__ << " key 0x" << std::hex << first_key << " bit 0x"
		 << i << " has 0x" << has << " expected 0x" << val << std::dec
		 << dendl;
	    ++errors;
	  }
	}
      } else {
	if (val) {
	  derr << __func__ << " key 0x" << std::hex << first_key
	       << " not present, expected 0x" << val << std::dec << dendl;
	  ++errors;
	}
      }
      first_key += bytes_per_key;
    }
    // middle keys
    if (first_key < last_key) {
      while (first_key < last_key) {
	string k;
	make_offset_key(first_key, &k);
	bufferlist bl;
	kvdb->get(bitmap_prefix, k, &bl);
	if (bl.length() > 0) {
	  const char *p = bl.c_str();
	  for (unsigned i = 0; i < blocks_per_key; ++i) {
	    int has = !!(p[i >> 3] & (1ull << (i & 7)));
	    if (has != val) {
	      derr << __func__ << " key 0x" << std::hex << first_key << " bit 0x"
		   << i << " has 0x" << has << " expected 0x" << val
		   << std::dec << dendl;
	      ++errors;
	    }
	  }
	} else {
	  if (val) {
	    derr << __func__ << " key 0x" << std::hex << first_key
		 << " not present, expected 0x" << val << std::dec << dendl;
	    ++errors;
	  }
	}
	first_key += bytes_per_key;
      }
    }
    ceph_assert(first_key == last_key);
    {
      string k;
      make_offset_key(first_key, &k);
      bufferlist bl;
      kvdb->get(bitmap_prefix, k, &bl);
      if (bl.length() > 0) {
	const char *p = bl.c_str();
	unsigned e = ((offset + length - 1) & ~key_mask) / bytes_per_block;
	for (unsigned i = 0; i < e; ++i) {
	  int has = !!(p[i >> 3] & (1ull << (i & 7)));
	  if (has != val) {
	    derr << __func__ << " key 0x" << std::hex << first_key << " bit 0x"
		 << i << " has 0x" << has << " expected 0x" << val << std::dec
		 << dendl;
	    ++errors;
	  }
	}
      } else {
	if (val) {
	  derr << __func__ << " key 0x" << std::hex << first_key
	       << " not present, expected 0x" << val << std::dec << dendl;
	  ++errors;
	}
      }
    }
  }
  if (errors) {
    derr << __func__ << " saw " << errors << " errors" << dendl;
    ceph_abort_msg("bitmap freelist errors");
  }
}

void BitmapFreelistManager::allocate(
  uint64_t offset, uint64_t length,
  KeyValueDB::Transaction txn)
{
  dout(10) << __func__ << " 0x" << std::hex << offset << "~" << length
	   << std::dec << dendl;
  _xor(offset, length, txn);
}

void BitmapFreelistManager::release(
  uint64_t offset, uint64_t length,
  KeyValueDB::Transaction txn)
{
  dout(10) << __func__ << " 0x" << std::hex << offset << "~" << length
	   << std::dec << dendl;
  _xor(offset, length, txn);
}

void BitmapFreelistManager::_xor(
  uint64_t offset, uint64_t length,
  KeyValueDB::Transaction txn)
{
  // must be block aligned
  ceph_assert((offset & block_mask) == offset);
  ceph_assert((length & block_mask) == length);

  uint64_t first_key = offset & key_mask;
  uint64_t last_key = (offset + length - 1) & key_mask;
  dout(20) << __func__ << " first_key 0x" << std::hex << first_key
	   << " last_key 0x" << last_key << std::dec << dendl;

  if (first_key == last_key) {
    bufferptr p(blocks_per_key >> 3);
    p.zero();
    unsigned s = (offset & ~key_mask) / bytes_per_block;
    unsigned e = ((offset + length - 1) & ~key_mask) / bytes_per_block;
    for (unsigned i = s; i <= e; ++i) {
      p[i >> 3] ^= 1ull << (i & 7);
    }
    string k;
    make_offset_key(first_key, &k);
    bufferlist bl;
    bl.append(p);
    dout(30) << __func__ << " 0x" << std::hex << first_key << std::dec << ": ";
    bl.hexdump(*_dout, false);
    *_dout << dendl;
    txn->merge(bitmap_prefix, k, bl);
  } else {
    // first key
    {
      bufferptr p(blocks_per_key >> 3);
      p.zero();
      unsigned s = (offset & ~key_mask) / bytes_per_block;
      unsigned e = blocks_per_key;
      for (unsigned i = s; i < e; ++i) {
	p[i >> 3] ^= 1ull << (i & 7);
      }
      string k;
      make_offset_key(first_key, &k);
      bufferlist bl;
      bl.append(p);
      dout(30) << __func__ << " 0x" << std::hex << first_key << std::dec << ": ";
      bl.hexdump(*_dout, false);
      *_dout << dendl;
      txn->merge(bitmap_prefix, k, bl);
      first_key += bytes_per_key;
    }
    // middle keys
    while (first_key < last_key) {
      string k;
      make_offset_key(first_key, &k);
      dout(30) << __func__ << " 0x" << std::hex << first_key << std::dec
      	 << ": ";
      all_set_bl.hexdump(*_dout, false);
      *_dout << dendl;
      txn->merge(bitmap_prefix, k, all_set_bl);
      first_key += bytes_per_key;
    }
    ceph_assert(first_key == last_key);
    {
      bufferptr p(blocks_per_key >> 3);
      p.zero();
      unsigned e = ((offset + length - 1) & ~key_mask) / bytes_per_block;
      for (unsigned i = 0; i <= e; ++i) {
	p[i >> 3] ^= 1ull << (i & 7);
      }
      string k;
      make_offset_key(first_key, &k);
      bufferlist bl;
      bl.append(p);
      dout(30) << __func__ << " 0x" << std::hex << first_key << std::dec << ": ";
      bl.hexdump(*_dout, false);
      *_dout << dendl;
      txn->merge(bitmap_prefix, k, bl);
    }
  }
}

//=====================================================================
//=====================================================================

/*******************************************************
 * function name: onebit_xor
 * function:save new bitmap as key-value into db
 * input: offset of pextent, count, txn
 * output: None
 * return: None
 * ****************************************************/
void BitmapFreelistManager::onebit_xor(uint64_t offset, KeyValueDB::Transaction txn, uint64_t min_alloc_size){
  std::lock_guard l(lock);
  dout(30) << __func__ << "====onebit_xor start. offset: " << offset << "=====" << dendl;
  //1.prepare parameter
  uint64_t chunk_size = min_alloc_size;
  uint64_t chunk_num = blocks_per_key;
  uint64_t base_offset = offset;

  //print and check
 // uint64_t chunk_index = offset / chunk_size / chunk_num;
 // uint64_t chunk_index_detail = ((offset / chunk_size) % chunk_num) / 8;
 // uint64_t chunk_index_detail_detail = ((offset / chunk_size) % chunk_num) % 8;
 // uint64_t index_offset = chunk_index_detail_detail;

 // unsigned int state_one = (0x01) & (bit_array_one[chunk_index][chunk_index_detail] >> index_offset);
 // unsigned int state_two = (0x01) & (bit_array_two[chunk_index][chunk_idnex_detail] >> index_offset);
 // dout(30) << __func__ << " before state_one: " << state_one << ", state_two: " << state_two << ", " << offset << dendl;
  //2.locate the index of target offset
  uint64_t base_number = base_offset / (chunk_num*chunk_size);
  dout(30) << __func__ << " base_offset: " << base_offset << ", base_number: " << base_number << dendl;

  //3.prepare key for index
  uint64_t key_offset = base_number * chunk_num * chunk_size;
  string k;
  k.reserve(10);
  _key_encode_u64_new(key_offset, &k);

  //4.prepare value for db
  bufferlist bl;
  uint64_t bytes_per_kv = chunk_num >> 2;
  bufferptr p(bytes_per_kv);
  p.zero();

  for(uint64_t i = 0; i < (bytes_per_kv/2); i++){
    p[i] = bit_array_one[base_number][i];
    p[i+16] = bit_array_two[base_number][i];
  }
  bl.append(p);

  //5.write in db
  txn->set(bitmap_prefix, k, bl);
  dout(30) << __func__ << "=====onebit_xor finish.=====" << dendl;
}

/**************************************************
 * function name: onebit_bitmap_create
 * function: create new bitmap during restart
 * input: pointer of kvdb
 * output: None
 * return: None
 * **********************************************/
int BitmapFreelistManager::onebit_bitmap_create(KeyValueDB *kvdb, uint64_t min_alloc_size){
  std::lock_guard l(lock);
  
  //0.pre alloc space
  uint64_t whole_block_size = size / min_alloc_size;
  uint64_t index_size = whole_block_size / blocks_per_key;

  bit_array_one.resize(index_size);
  bit_array_two.resize(index_size);
  for(uint64_t i = 0; i < index_size; i++){
    bit_array_one[i] = "0000000000000000";
    bit_array_two[i] = "0000000000000000";
    for(int j = 0; j < 16; j++){
      bit_array_one[i][j] &= 0x0;
      bit_array_two[i][j] &= 0x0;
    }
  }
  dout(30) << __func__ << " whole_block_size: " << whole_block_size << ", index_size: " << index_size << dendl;

  //1.firstly, get the pointer of kvdb
  dout(30) << __func__ << " ====let's create onebit_bitmap for counting reuse pextents.====" << dendl;
  enumerate_p = kvdb->get_iterator(bitmap_prefix);
  enumerate_p->lower_bound(string());

  //2.create new bitmap
  while(enumerate_p->valid()){
    string k = enumerate_p->key();
    dout(30) << "====key size: " << k.size() << ", k: " << k << dendl;

    if(k.length() == 10){
      //2.1 read and convert key
      const char *p = k.c_str();
      uint64_t my_offset;
      p = _key_decode_u64_new(p, &my_offset);
      dout(30) << "====1.my_offset: " << my_offset << dendl;

      //2.2 then get value according to key, and build as count_array in memory
      enumerate_bl = enumerate_p->value();
      const char* segment_p = enumerate_bl.c_str();
      dout(30) << "enumerate_bl: " << enumerate_bl << ", enumerate_bl length: " << enumerate_bl.length() << dendl;
      string segment_value = "00000000000000000000000000000000";
      for(int i = 0; i < 32; i++){
        segment_value[i] = segment_p[i];
      }
      dout(30) << "===2.segment_value:" << segment_value << ", segment_value length:" << segment_value.length() << dendl;

      uint64_t array_pos = my_offset / min_alloc_size / blocks_per_key;
      bit_array_one[array_pos] = segment_value.substr(0,16);
      bit_array_two[array_pos] = segment_value.substr(16);
      dout(30) << __func__ << " my_offset:" << my_offset << ", array_pos:" << array_pos
               << ", bit_array_one:" << bit_array_one[array_pos]
               << ", bit_array_two:" << bit_array_two[array_pos] << dendl;
    }

    //3.prepare for next kv
    enumerate_p->next();
    if(!enumerate_p->valid())
      break;
  }
  dout(30) << __func__ << "====onebit_bitmap has already created.====" << dendl;
  dout(30) << __func__ << "reuse_count: " << reuse_count() << dendl;
  return 0;
}

/***********************************************************
 * fucntion name: onebit_check_bitmap
 * function: check during restart
 * input: extents, compressed, min_alloc_size
 * output: None
 * return: None
 * *********************************************************/
void BitmapFreelistManager::onebit_check_bitmap(const PExtentVector& extents, bool compressed, uint64_t min_alloc_size){
  if(compressed){
    dout(30) << __func__ << " extents.size: " << extents.size() << dendl;

    for(auto e:extents){
      dout(30) << __func__ << " offset: " << e.offset << " length: " << e.length << dendl;

      //1.first align offset and length
     uint64_t chunk_size = min_alloc_size;
     uint64_t new_start = p2align(e.offset, chunk_size);
     uint64_t new_end = p2align(e.offset + e.length-1, chunk_size);
     dout(30) << __func__ << "After p2align, offset: " << new_start << ", end(offset+length): " << new_end << dendl;

     //2.check after align
     if(new_start == new_end){
       uint64_t pos_index = new_start / chunk_size;
       if(pos_index < count_number.size()){
         dout(30) << __func__ << "1.num before check: " << int(count_number[pos_index]) << dendl;
         count_number[pos_index]--;
         if(count_number[pos_index] < 0){
           dout(30) << "ERROR!" << dendl;
         }
       }
       else{
         dout(30) << __func__ << "offset is not reused." << dendl;
         return;
       }
     }
     else{
       uint64_t temp_offset = new_start;
       while(temp_offset <= new_end){
         uint64_t pos_index = temp_offset / chunk_size;
         if(pos_index < count_number.size()){
           dout(30) << __func__ << "2.num before check: " << int(count_number[pos_index]) << ", offset: " << temp_offset << dendl;
           count_number[pos_index]--;
           if(count_number[pos_index] < 0){
             dout(30) << "ERROR!" << dendl;
           }
           temp_offset += chunk_size;
         }
         else{
           break;
         }
       }
     }
    }
  }
}

/*********************************************************
 * function name: convert_num
 * function: from bit to num
 * input: None
 * output: None
 * return: None
 * ******************************************************/
void BitmapFreelistManager::convert_num(){
  uint64_t index_count = bit_array_one.size();
  for (uint64_t i = 0; i < index_count; i++) {
    for (int j = 0; j < 16; j++) {
      for(unsigned int k = 0; k < 8; k++) {
        unsigned int state_one = (0x1) & (bit_array_one[i][j] >> k);
        unsigned int state_two = (0x1) & (bit_array_two[i][j] >> k);
        int8_t result = static_cast<int8_t>((state_one != 0)) + 2 * static_cast<uint8_t>((state_two != 0));
        count_number.push_back(result);
      }
    }
  }
  dout(30) << __func__ << " count array has already be created. number: " << count_number.size() << dendl;
}

/********************************************************
 * function name: reuse_count
 * function: check during restart
 * input: None
 * output: None
 * return: None
 * *****************************************************/
int64_t BitmapFreelistManager::reuse_count(){
  return ReuseCountClose(count_number);
}

/*******************************************************
 * function name: plus_function
 * function: count++
 * input: offset
 * output: None
 * return None
 * ***************************************************/
void BitmapFreelistManager::plus_function(uint64_t offset, uint64_t min_alloc_size){
  std::lock_guard l(lock);
  if (min_alloc_size == 0 || blocks_per_key == 0) {
    derr << __func__ << "unrecognized min_alloc_size or blocks_per_key " << dendl;
  }

  uint64_t chunk_index = offset / min_alloc_size /blocks_per_key;
  uint64_t chunk_index_detail = ((offset / min_alloc_size) % blocks_per_key) / 8;
  uint64_t chunk_index_detail_internal = ((offset / min_alloc_size) % blocks_per_key) % 8;
  uint64_t index_offset = chunk_index_detail_internal;

  unsigned int state_one = (0x01) & (bit_array_one[chunk_index][chunk_index_detail] >> index_offset);
  unsigned int state_two = (0x01) & (bit_array_two[chunk_index][chunk_index_detail] >> index_offset);

  if ((state_one == 0) && (state_two == 0)) {
    bit_array_one[chunk_index][chunk_index_detail] |= (1 << index_offset);
  } else if ((state_one != 0) && (state_two == 0)) {
    bit_array_one[chunk_index][chunk_index_detail] &= ~(1 << index_offset);
    bit_array_two[chunk_index][chunk_index_detail] |= (1 << index_offset);
  } else if ((state_one == 0) && (state_two != 0)) {
    bit_array_one[chunk_index][chunk_index_detail] |= (1 << index_offset);
  } else {
    derr << __func__ << " something wrong happened: state_one or state_two in wrong station " << dendl;
  }

  state_one = (0x01) & (bit_array_one[chunk_index][chunk_index_detail] >> index_offset);
  state_two = (0x01) & (bit_array_two[chunk_index][chunk_index_detail] >> index_offset);
}

/**************************************************************
 * function name: minus_function
 * function: count--
 * input: offset
 * output: None
 * return: None
 * **********************************************************/
void BitmapFreelistManager::minus_function(uint64_t offset, uint64_t min_alloc_size){
  std::lock_guard l(lock);
  if (min_alloc_size == 0 || blocks_per_key == 0) {
    derr << __func__ << " unrecognized min_alloc_size or blocks_per_key " << dendl; 
  }

  uint64_t chunk_index = offset / min_alloc_size /blocks_per_key;
  uint64_t chunk_index_detail = ((offset / min_alloc_size) % blocks_per_key) / 8;
  uint64_t chunk_index_detail_internal = ((offset / min_alloc_size) % blocks_per_key) % 8;
  uint64_t index_offset = chunk_index_detail_internal;

  unsigned int state_one = (0x01) & (bit_array_one[chunk_index][chunk_index_detail] >> index_offset);
  unsigned int state_two = (0x01) & (bit_array_two[chunk_index][chunk_index_detail] >> index_offset);

  if ((state_one != 0 ) && (state_two == 0)) {
    bit_array_one[chunk_index][chunk_index_detail] &= ~(1 << index_offset);
  } else if ((state_one == 0) && (state_two != 0)) {
    bit_array_one[chunk_index][chunk_index_detail] |= (1 << index_offset);
    bit_array_two[chunk_index][chunk_index_detail] &= ~(1 << index_offset);
  } else if ((state_one != 0) && (state_two != 0)) {
    bit_array_one[chunk_index][chunk_index_detail] &= ~(1 << index_offset);
  } else {
    derr << __func__ << " something wrong happened: state_one or state_two in wrong station " << dendl;
  }

  state_one = (0x01) & (bit_array_one[chunk_index][chunk_index_detail] >> index_offset);
  state_two = (0x01) & (bit_array_two[chunk_index][chunk_index_detail] >> index_offset);
}

/*******************************************************
 * function name: zero_check_function
 * function: count check
 * input: offset, min_alloc_size
 * output: None
 * return: 0 if count=0, 1 if count !=0
 * ***************************************************/
int BitmapFreelistManager::zero_check_function(uint64_t offset, uint64_t min_alloc_size){
  if (min_alloc_size == 0 || blocks_per_key == 0) {
    derr << __func__ << " unrecognized min_alloc_size or blocks_per_key " << dendl;  
  }
  uint64_t chunk_index = offset / min_alloc_size / blocks_per_key;
  uint64_t chunk_index_detail = ((offset / min_alloc_size) % blocks_per_key) / 8;
  uint64_t chunk_index_detail_internal = ((offset / min_alloc_size) % blocks_per_key) % 8;
  uint64_t index_offset = chunk_index_detail_internal;

  unsigned int state_one = (0x01) & (bit_array_one[chunk_index][chunk_index_detail] >> index_offset);
  unsigned int state_two = (0x01) & (bit_array_two[chunk_index][chunk_index_detail] >> index_offset);

  if ((state_one == 0)&& (state_two == 0)) {
    return 0;
  } else {
    return 1;
  }
}
