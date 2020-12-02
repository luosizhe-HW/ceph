#ifndef HI_CORE_UTIL_H
#define HI_CORE_UTIL_H

#include <stdint.h>
#include <map>
#include <vector>
#include <set>
#include <boost/tuple/tuple.hpp>
using namespace std;

bool HiSetWriteSection( const uint64_t offest, const uint64_t len, const uint64_t csize,
		uint64_t &headStart, uint64_t &headLength);

bool HiRebuildToread( const map<uint64_t, uint64_t> &writeSet, uint64_t csize,
		map<uint64_t, uint64_t> &toRead);

struct HiECInfo {
	size_t dataCnt;
	size_t chunkCnt;
	const vector<int> &cMapping;
	uint64_t cSize;
	uint64_t sWidth;

	HiECInfo(size_t k,size_t km, const vector<int> &chunkMapping,uint64_t chunkSize,
			uint64_t stripeWidth)
		: dataCnt(k),
		  chunkCnt(km),
		  cMapping(chunkMapping),
		  cSize(chunkSize),
		  sWidth(stripeWidth) {}
};

void HiGetRelatedShards(const pair<uint64_t,uint64_t> &rRange,
			const HiECInfo &ecInfo, set<int> &out);

void HiGetShardsRangeeToRead(const pair<uint64_t,uint64_t> &rRange,
			const HiECInfo &ecInfo,
			map<int , set<pair<uint64_t, uint64_t>>> &wants);

void HiGetWriteToshards(const map<uint64_t,uint64_t> &toWriteChunkAlign,
			const HiECInfo &ecInfo ,set<int> &wantToWrite);

void HiGetReconstructShards(const unsigned int start,const unsigned int count,
			const uint64_t len, const HiECInfo &ecInfo,
			vector<int> &shardID);

void HiReconstructPrepare(
		const HiECInfo &ecInfo, const pair<uint64_t,uint64_t> &rRange,
		vector<boot::tuple<unsigned int , unsigned int , unsigned int>> &cInfo);

#endif











