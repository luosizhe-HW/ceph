// -*- mode:C++; tab-width:8; c-basic-offset:2; indent-tabs-mode:t -*-
// vim: ts=8 sw=2 smarttab

#ifndef CEPH_OS_BLUESTORE_FREELISTMANAGER_H
#define CEPH_OS_BLUESTORE_FREELISTMANAGER_H

#include <string>
#include <map>
#include <mutex>
#include <ostream>
#include "kv/KeyValueDB.h"

#include "BlueStore.h"

class FreelistManager {
public:
  map<uint64_t, int8_t> new_bitmap;
  vector<string> count_array;
  vector<string> bit_array_one;
  vector<string> bit_array_two;
  vector<int8_t> count_number;

  CephContext* cct;
  FreelistManager(CephContext* cct) : cct(cct) {}
  virtual ~FreelistManager() {}

  static FreelistManager *create(
    CephContext* cct,
    string type,
    string prefix);

  static void setup_merge_operators(KeyValueDB *db);

  virtual int create(uint64_t size, uint64_t granularity,
		     KeyValueDB::Transaction txn) = 0;

  virtual int expand(uint64_t new_size,
		     KeyValueDB::Transaction txn) = 0;

  virtual int init(KeyValueDB *kvdb) = 0;
  virtual void shutdown() = 0;

  virtual void dump(KeyValueDB *kvdb) = 0;

  virtual void enumerate_reset() = 0;
  virtual bool enumerate_next(KeyValueDB *kvdb, uint64_t *offset, uint64_t *length) = 0;

  virtual void allocate(
    uint64_t offset, uint64_t length,
    KeyValueDB::Transaction txn) = 0;
  virtual void release(
    uint64_t offset, uint64_t length,
    KeyValueDB::Transaction txn) = 0;

  virtual uint64_t get_size() const = 0;
  virtual uint64_t get_alloc_units() const = 0;
  virtual uint64_t get_alloc_size() const = 0;

  virtual void onebit_xor(uint64_t offset, KeyValueDB::Transaction txn, uint64_t min_alloc_size) = 0;
  virtual int onebit_bitmap_create(KeyValueDB *kvdb, uint64_t min_alloc_size) = 0;
  virtual void onebit_check_bitmap(const PExtentVector& extents, bool compressed, uint64_t min_alloc_size) = 0;
  virtual void convert_num() = 0;
  virtual int64_t reuse_count() = 0;
  virtual void plus_function(uint64_t offset, uint64_t min_alloc_size) = 0;
  virtual void minus_function(uint64_t offset, uint64_t min_alloc_size) = 0;
  virtual int zero_check_function(uint64_t offset, uint64_t min_alloc_size) = 0;
};

#endif
