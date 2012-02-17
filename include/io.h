/*
 * io.h
 */
#ifndef IO_H
#define IO_H

#include "types.hpp"

namespace io {

bool readFile(types::String filename, types::String& content);

bool readBinary(types::String filename, types::Binary& content);
bool writeBinary(types::String filename, const types::Binary &content);
bool writeStream(types::String filename, const types::Stream &content);

}
#endif
