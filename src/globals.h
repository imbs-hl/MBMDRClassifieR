
#ifndef GLOBALS_H_
#define GLOBALS_H_

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
TypeName(const TypeName&);                 \
void operator=(const TypeName&)

typedef unsigned int uint;

// Default values
const uint DEFAULT_NUM_THREADS = 0;

#endif /* GLOBALS_H_ */
