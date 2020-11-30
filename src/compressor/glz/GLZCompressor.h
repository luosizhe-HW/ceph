// -*- mode:C++; tab-width:8; c-basic-offset:2; indent-tabs-mode:t -*-
// vim: ts=8 sw=2 smarttab

#ifndef CEPH_GLZCOMPRESSOR_H
#define CEPH_GLZCOMPRESSOR_H

#include <glz.h>

#include "compressor/Compressor.h"
#include "include/buffer.h"
#include "include/encoding.h"
#include "common/config.h"
#include "common/Tub.h"

#define GLZ_DEFAULT_LEVEL 1

class GLZCompressor : public Compressor {
 CephContext *const cct;
 public:
  GLZCompressor(CephContext* cct) : Compressor(COMP_ALG_GLZ, "glz"), cct(cct) {}
  ~GLZCompressor() {}

  int compress(const bufferlist &src, bufferlist &dst, boost::optional<int32_t> &compressor_message) override {
    if (!src.is_contiguous()) {
      bufferlist new_src = src;
      new_src.rebuild();
      return compress(new_src, dst, compressor_message);
    }
    bufferptr outptr = buffer::create_small_page_aligned(src.length());

    auto p = src.begin();
    size_t left = src.length();
    int pos = 0;
    const char *data = nullptr;
    unsigned num = src.get_num_buffers();
    encode((uint32_t)num, dst);
    glz_encoder* glz_stream = nullptr;
    glz_stream = glz_compress_initial((long unsigned)left);
    if (glz_stream == nullptr) {
	    return -2;
    }
    encode_model glz_mode = GLZ_FAST;
    if (cct->_conf->compressor_glz_level == GLZ_HIGH_CR) {
	    glz_mode = GLZ_HIGH_CR;
    }
    while (left) {
      uint32_t origin_len = p.get_ptr_and_advance(left, &data);
      int compressed_len = glz_compress_default((glz_encoder *const)glz_stream,outptr.c_str()+pos,outptr.length()-pos,data,origin_len,glz_mode,GLZ_DEFAULT_LEVEL);
      if (compressed_len <= 0) {
	glz_compress_delete(&glz_stream);
        return -1;
    }
      pos += compressed_len;
      left -= origin_len;
      encode(origin_len, dst);
      encode((uint32_t)compressed_len, dst);
  }
  ceph_assert(p.end());
  compressor_message = cct->_conf->compressor_glz_level;
  dst.append(outptr, 0, pos);
  glz_compress_delete(&glz_stream);
 return 0;
} 

  int decompress(const bufferlist &src, bufferlist &dst, boost::optional<int32_t> compressor_message) override {
    auto i = std::cbegin(src);
    return decompress(i, src.length(), dst, compressor_message);
  }

  int decompress(bufferlist::const_iterator &p,
		 size_t compressed_len,
		 bufferlist &dst,
		 boost::optional<int32_t> compressor_message) override {
	  uint32_t count;
	  std::vector<std::pair<uint32_t,uint32_t> > compressed_pairs;
	  decode(count, p);
	  compressed_pairs.resize(count);
	  uint32_t total_origin = 0;
	  for (unsigned i = 0;i < count; ++i) {
	  	decode(compressed_pairs[i].first, p);
		decode(compressed_pairs[i].second, p);
		total_origin += compressed_pairs[i].first;
	  }
	  compressed_len -= (sizeof(uint32_t) + sizeof(uint32_t) * count * 2);
	  bufferptr dstptr(total_origin);
	  bufferptr cur_ptr = p.get_current_ptr();
	  bufferptr *ptr = &cur_ptr;
	  Tub<bufferptr> data_holder;
	  if (compressed_len != cur_ptr.length()) {
	    data_holder.construct(compressed_len);
	    p.copy_deep(compressed_len, *data_holder);
	    ptr = data_holder.get();
	  }
	  char *c_in = ptr->c_str();
	  char *c_out = dstptr.c_str();
	  encode_model glz_mode = GLZ_FAST;
	  if (*compressor_message == GLZ_HIGH_CR) {
	    glz_mode = GLZ_HIGH_CR;
	  }
	  for (unsigned i = 0; i < count; ++i) {
		  int r = glz_decompress_default(c_in, c_out, compressed_pairs[i].first, c_in, compressed_pairs[i].second, glz_mode);
		  if (r == (int)compressed_pairs[i].first) {
			  c_in += compressed_pairs[i].second;
			  c_out += compressed_pairs[i].first;
		  } else if (r < 0) {
			  return -1;
		  } else { 
			  return -2;
		  }
	  }
	  dst.push_back(std::move(dstptr));
	  return 0;
  }
};

#endif
