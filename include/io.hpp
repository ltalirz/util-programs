/**
 * @file io.hpp
 */
#ifndef IO_HPP
#define IO_HPP

#include "types.hpp"

namespace io {

bool readFile(types::String filename, types::String& content);

bool readBinary(types::String filename, types::Binary& content);
bool writeBinary(types::String filename, const types::Binary &content);
bool writeStream(types::String filename, const types::Stream &content);

types::String getAbsolutePath(types::String relativePath);
types::String getFileName(types::String path);
types::String getFileExtension(types::String path);

}
#endif
