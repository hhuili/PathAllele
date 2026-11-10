# WFA2 interface headers

Unmodified include closure of `wavefront/wfa.hpp` from pjotrp/WFA2-lib commit
`f86a35abf89933605d81df6e67e42edc0aaf7f55`, matching the previously validated backend.
Source: https://github.com/pjotrp/WFA2-lib/tree/f86a35abf89933605d81df6e67e42edc0aaf7f55
Archive SHA256: `901070d8cc5ebaafb372a0b1bb95edd164a10e3db4e50d1da15d9f3116e8e73a`.

Only interface headers and their upstream MIT license are included. The alignment
implementation is supplied by the linked vcflib/WFA2 runtime package. These headers
are needed because the tested Linux vcflib package ships the libraries but not the
complete WFA2 include tree. Do not update independently of backend ABI/parity tests.
