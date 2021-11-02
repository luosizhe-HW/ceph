#pragma once

#include <mutex>
#include <boost/intrusive/avl_set.hpp>

#include "Allocator.h"
#include "os/bluestore/bluestore_types.h"
#include "include/mempool.h"

struct range_seg_t {
  MEMPOOL_CLASS_HELPERS();
  uint64_t start;
  uint64_t end;

  range_seg_t(uint64_t start, uint64_t end)
    : start{start},
      end{end}
    {}

  struct before_t {
    template<typename KeyLeft, typename KeyRight>
    bool operator()(const KeyLeft& lhs, const KeyRight& rhs) const {
      return lhs.end <= rhs.start;
    }
  };
  boost::intrusive::avl_set_member_hook<> offset_hook;

  struct shorter_t {
    template<typename KeyType>
    bool operator()(const range_seg_t& lhs, const KeyType& rhs) const {
      auto lhs_size = lhs.end - lhs.start;
      auto rhs_size = rhs.end - rhs.start;
      if (lhs_size < rhs_size) {
          return true;
      } else if (lhs_size > rhs_size) {
          return false;
      } else {
          return lhs.start < rhs.start;
      }
    }
  };
  inline uint64_t length() const {
    return end - start;
  }
  boost::intrusive::avl_set_member_hook<> size_hook;
};

class AvlAllocator : public Allocator {
  struct dispose_rs {
    void operator()(range_seg_t* p){
      delete p;
    }
  };

protected:
  AvlAllocator(CephContext* cct, int64_t device_size, int64_t block_size,
    uint64_t max_mem, 
    const std::string& name);

public:
  AvlAllocator(CephContext* cct, int64_t device_size, int64_t block_size,
    const std::string& name);
  ~AvlAllocator();
  int64_t allocate(
    uint64_t want,
    uint64_t unit,
    uint64_t max_alloc_size,
    int64_t hint,
    PExtentVector *extents) override;
  void release(const interval_set<uint64_t>& release_set) override;
  int64_t get_capacity() const {
    return num_total;
  }

  uint64_t get_block_size() const {
    return block_size;
  }
  uint64_t get_free() override;
  double get_fragmentation();

  void dump() override;
  void dump(std::function<void(uint64_t offset, uint64_t length)> notify) override;
  void init_add_free(uint64_t offset, uint64_t length) override;
  void init_rm_free(uint64_t offset, uint64_t length) override;
  void shutdown()  override;

private:
  template<class Tree>
  uint64_t _block_picker(const Tree& t, uint64_t *cursor, uint64_t size,
    uint64_t align);
  int _allocate(
    uint64_t size,
    uint64_t unit,
    uint64_t *offset,
    uint64_t *length);

  using range_tree_t = 
    boost::intrusive::avl_set<
      range_seg_t,
      boost::intrusive::compare<range_seg_t::before_t>,
      boost::intrusive::member_hook<
    	range_seg_t,
    	boost::intrusive::avl_set_member_hook<>,
    	&range_seg_t::offset_hook>>;
  range_tree_t range_tree;    ///< main range tree

  using range_size_tree_t=
    boost::intrusive::avl_multiset<
      range_seg_t,
      boost::intrusive::compare<range_seg_t::shorter_t>,
      boost::intrusive::member_hook<
  range_seg_t,
  boost::intrusive::avl_set_member_hook<>,
  &range_seg_t::size_hook>,
      boost::intrusive::constant_time_size<true>>;
  range_size_tree_t range_size_tree;
    
  const int64_t num_total;   ///<device size
  const uint64_t block_size; ///<block size
  uint64_t num_free = 0;     ///< total bytes in freelist

  static constexpr unsigned MAX_LBAS = 64;
  uint64_t lbas[MAX_LBAS] = {0};
  uint64_t range_size_alloc_threshold = 0;
  int range_size_alloc_free_pct = 0;
  uint64_t range_count_cap = 0;

  void _range_size_tree_rm(range_seg_t& r) {
    ceph_assert(num_free >= r.length());
    num_free -= r.length();
    range_size_tree.erase(r);
  }
  void _range_size_tree_try_insert(range_seg_t& r) {
    if (_try_insert_range(r.start, r.end)) {
      range_size_tree.insert(r);
      num_free += r.length();
    } else {
      range_tree.erase_and_dispose(r, dispose_rs{});
    }
  }
  bool _try_insert_range(uint64_t start, 
			 uint64_t end,
			 range_tree_t::iterator* insert_pos = nullptr) {
    bool res = !range_count_cap || range_size_tree.size() < range_count_cap;
    bool remove_lowest = false;
    if (!res) {
      if (end - start > _lowest_size_available()) {
        remove_lowest = true;
	res = true;
      }
    }
    if (!res) {
      _spillover_range(start, end);
    } else {
      if (insert_pos) {
        auto new_rs = new range_seg_t{start, end};
	range_tree.insert_before(*insert_pos, *new_rs);
	range_size_tree.insert(*new_rs);
	num_free += new_rs->length();
      }
      if (remove_lowest) {
        auto r = range_size_tree.begin();
	_range_size_tree_rm(*r);
	_spillover_range(r->start, r->end);
	range_tree.erase_and_dispose(*r, dispose_rs{});
      }
    }
    return res;
  }
  virtual void _spillover_range(uint64_t start, uint64_t end) {
    ceph_assert(false);
  }

protected:
  virtual void _add_to_tree(uint64_t start, uint64_t size);

protected:
  CephContext* cct;
  std::mutex lock;

  double _get_fragmentation() const {
    auto free_blocks = p2align(num_free, block_size) / block_size;
    if (free_blocks <= 1) {
      return .0;
    }
    return (static_cast<double>(range_tree.size() - 1) / (free_blocks - 1));
  }
  void _dump() const;

  uint64_t _lowest_size_available() {
    auto rs = range_size_tree.begin();
    return rs != range_size_tree.end() ? rs->length() : 0;
  }

  int64_t _allocate(
    uint64_t want,
    uint64_t unit,
    uint64_t max_alloc_size,
    int64_t hint,
    PExtentVector *extents);

  void _release(const interval_set<uint64_t>& release_set);
  void _release(const PExtentVector& release_set);
  void _shutdown();

  void _process_range_removal(uint64_t start, uint64_t end, range_tree_t::iterator& rs);
  void _remove_from_tree(uint64_t start, uint64_t size);
  void _try_remove_from_tree(uint64_t start, uint64_t size,
    std::function<void(uint64_t offset, uint64_t length, bool found)> cb);

  uint64_t _get_free() const {
    return num_free;
  }
};

