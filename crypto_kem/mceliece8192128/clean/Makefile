# This Makefile can be used with GNU Make or BSD Make

<<<<<<< HEAD
LIB = libmceliece8192128_clean.a

SOURCES = aes256ctr.c benes.c bm.c controlbits.c crypto_int16.c crypto_int32.c \
		crypto_uint16.c crypto_uint32.c crypto_uint64.c decrypt.c encrypt.c gf.c \
		operations.c pk_gen.c root.c sk_gen.c synd.c transpose.c util.c

HEADERS = aes256ctr.h api.h api.h benes.h bm.h controlbits.h crypto_hash.h crypto_hash.h \
=======
# Defined variables with reusable file names, compiler flags, and targets

# The library to be created.
LIB = libmceliece8192128_clean.a

SOURCES = ../../../common/randombytes.c ../../../common/fips202.c ../../../common/sha2.c aes256ctr.c benes.c bm.c controlbits.c crypto_int16.c crypto_int32.c \
		crypto_uint16.c crypto_uint32.c crypto_uint64.c decrypt.c encrypt.c gf.c \
		operations.c pk_gen.c root.c sk_gen.c synd.c transpose.c util.c

HEADERS = ../../../common/randombytes.h ../../../common/fips202.h ../../../common/sha2.h aes256ctr.h api.h api.h benes.h bm.h controlbits.h crypto_hash.h crypto_hash.h \
>>>>>>> master
		crypto_int16.h crypto_int32.h crypto_kem.h crypto_uint16.h crypto_uint32.h \
		crypto_uint64.h decrypt.h encrypt.h gf.h int32_sort.h int32_sort.h namespace.h \
		operations.h params.h pk_gen.h root.h sk_gen.h synd.h transpose.h uint64_sort.h \
		util.h

<<<<<<< HEAD
OBJECTS = aes256ctr.o benes.o bm.o controlbits.o crypto_int16.o crypto_int32.o \
		crypto_uint16.o crypto_uint32.o crypto_uint64.o decrypt.o encrypt.o gf.o \
		operations.o pk_gen.o root.o sk_gen.o synd.o transpose.o util.o
=======
OBJECTS = ../../../common/randombytes.o ../../../common/fips202.o ../../../common/sha2.o aes256ctr.o benes.o bm.o controlbits.o crypto_int16.o crypto_int32.o \
          crypto_uint16.o crypto_uint32.o crypto_uint64.o decrypt.o encrypt.o gf.o \
          operations.o pk_gen.o root.o sk_gen.o synd.o transpose.o util.o
>>>>>>> master

CFLAGS  = -O3 -std=c99 -Wall -Wextra -pedantic -Werror -Wpedantic \
	  -Wredundant-decls -Wcast-align -Wmissing-prototypes \
	  -I../../../common/ $(EXTRAFLAGS)

<<<<<<< HEAD
all: $(LIB)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS)   -c -o $@ $<

$(LIB): $(OBJECTS)
	$(AR) -r $@ $(OBJECTS)

=======
# Rule builds the library. 
all: $(LIB)

# Pattern rule to compile .c files into .o files. $@ represents the target file (aes256ctr.o)
# $< represents the first prerequisite (aes256ctr.c). The $(CC) is the compiler, and

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS)   -c -o $@ $<

# LIB is the target  to be created. OBJECTS are the dependencies needed to create the 
# target. $(AR) -r $@ $(OBJECTS) commands are the steps to create the target.
# Archive the .o files into a static library.
$(LIB): $(OBJECTS)
	$(AR) -r $@ $(OBJECTS)

# Commast for specific tasks. Clean removed the specific filed and the library.
>>>>>>> master
clean:
	$(RM) $(OBJECTS)
	$(RM) $(LIB)