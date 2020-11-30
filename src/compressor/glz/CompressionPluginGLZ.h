#ifndef CEPH_COMPRESSION_PLUGIN_GLZ_H
#define CEPH_COMPRESSION_PLUGIN_GLZ_H

// -----------------------------------------------------------------------------
#include "ceph_ver.h"
#include "compressor/CompressionPlugin.h"
#include "GLZCompressor.h"
// -----------------------------------------------------------------------------

class CompressionPluginGLZ : public CompressionPlugin {

public:

  explicit CompressionPluginGLZ(CephContext* cct) : CompressionPlugin(cct)
  {}
  ~CompressionPluginGLZ() {}

  int factory(CompressorRef *cs, std::ostream *ss) override {
    if (compressor == 0) {
      GLZCompressor *interface = new GLZCompressor(cct);
      compressor = CompressorRef(interface);
    }
    *cs = compressor;
    return 0;
  }
};

#endif
