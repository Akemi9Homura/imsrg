"""Dependency-free constants shared by MR reference readers and writers."""

MAGIC = b"mrimsrg_jref1\0\0\0"
ENDIAN_MARKER = 0x01020304
SCHEMA_VERSION = 1
# Version 2 appends the optional reference energy and lambda3 blocks; files
# without either are still written as version 1.
SCHEMA_VERSION_LAMBDA3 = 2
SUPPORTED_SCHEMA_VERSIONS = (1, 2)
REFERENCE_ENERGY_LABEL_BYTES = 64
