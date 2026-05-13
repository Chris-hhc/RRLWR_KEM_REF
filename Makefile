CC = gcc
CFLAGS += -O3
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -fomit-frame-pointer
RM = /bin/rm

BASE_FOLDER = .
ARITH_FOLDER = ./arith
UTILS_FOLDER = ./utils
SOURCES_RRLWR = $(ARITH_FOLDER)/fprime.c $(ARITH_FOLDER)/poly.c $(ARITH_FOLDER)/ring.c $(ARITH_FOLDER)/packing.c $(ARITH_FOLDER)/uniform.c pke.c kem.c
SOURCES = $(SOURCES_RRLWR) $(UTILS_FOLDER)/auxfunc.c
HEADERS_RRLWR = $(ARITH_FOLDER)/fprime.h $(ARITH_FOLDER)/poly.h $(ARITH_FOLDER)/ring.h $(ARITH_FOLDER)/packing.h $(ARITH_FOLDER)/uniform.h parameters.h pke.h kem.h
HEADERS = $(HEADERS_RRLWR) $(UTILS_FOLDER)/auxfunc.h

.PHONY: all test clean

all: test speed KAT

test: \
  test/test_KEM128 \
  test/test_KEM256 \
  test/test_KEM512 \
  test/unit_tests_KEM128 \
  test/unit_tests_KEM256 \
  test/unit_tests_KEM512 \

speed: \
  test/test_speed_KEM128 \
  test/test_speed_KEM256 \
  test/test_speed_KEM512 \

KAT: \
  KAT/KAT_KEM128 \
  KAT/KAT_KEM256 \
  KAT/KAT_KEM512 \

test/test_KEM128: $(SOURCES) $(HEADERS) $(UTILS_FOLDER)/drng.c $(UTILS_FOLDER)/drng.h test/test_KEM.c
	$(CC) -I$(BASE_FOLDER) -I$(UTILS_FOLDER) -I$(ARITH_FOLDER) -DRRLWR_SECURITY_LEVEL=128 $(CFLAGS) $(SOURCES) $(UTILS_FOLDER)/drng.c test/test_KEM.c -o $@

test/test_KEM256: $(SOURCES) $(HEADERS) $(UTILS_FOLDER)/drng.c $(UTILS_FOLDER)/drng.h test/test_KEM.c
	$(CC) -I$(BASE_FOLDER) -I$(UTILS_FOLDER) -I$(ARITH_FOLDER) -DRRLWR_SECURITY_LEVEL=256 $(CFLAGS) $(SOURCES) $(UTILS_FOLDER)/drng.c test/test_KEM.c -o $@

test/test_KEM512: $(SOURCES) $(HEADERS) $(UTILS_FOLDER)/drng.c $(UTILS_FOLDER)/drng.h test/test_KEM.c
	$(CC) -I$(BASE_FOLDER) -I$(UTILS_FOLDER) -I$(ARITH_FOLDER) -DRRLWR_SECURITY_LEVEL=512 $(CFLAGS) $(SOURCES) $(UTILS_FOLDER)/drng.c test/test_KEM.c -o $@

test/unit_tests_KEM128: $(SOURCES) $(HEADERS) $(UTILS_FOLDER)/drng.c $(UTILS_FOLDER)/drng.h test/unit_tests_KEM.c
	$(CC) -I$(BASE_FOLDER) -I$(UTILS_FOLDER) -I$(ARITH_FOLDER) -DRRLWR_SECURITY_LEVEL=128 $(CFLAGS) $(SOURCES) $(UTILS_FOLDER)/drng.c test/unit_tests_KEM.c -o $@

test/unit_tests_KEM256: $(SOURCES) $(HEADERS) $(UTILS_FOLDER)/drng.c $(UTILS_FOLDER)/drng.h test/unit_tests_KEM.c
	$(CC) -I$(BASE_FOLDER) -I$(UTILS_FOLDER) -I$(ARITH_FOLDER) -DRRLWR_SECURITY_LEVEL=256 $(CFLAGS) $(SOURCES) $(UTILS_FOLDER)/drng.c test/unit_tests_KEM.c -o $@

test/unit_tests_KEM512: $(SOURCES) $(HEADERS) $(UTILS_FOLDER)/drng.c $(UTILS_FOLDER)/drng.h test/unit_tests_KEM.c
	$(CC) -I$(BASE_FOLDER) -I$(UTILS_FOLDER) -I$(ARITH_FOLDER) -DRRLWR_SECURITY_LEVEL=512 $(CFLAGS) $(SOURCES) $(UTILS_FOLDER)/drng.c test/unit_tests_KEM.c -o $@

test/test_speed_KEM128: $(SOURCES) $(HEADERS) test/cpucycles.h test/cpucycles.c test/speed_print.c test/speed_print.h $(UTILS_FOLDER)/drng.c $(UTILS_FOLDER)/drng.h test/test_speed_KEM.c
	$(CC) -I$(BASE_FOLDER) -I$(UTILS_FOLDER) -I$(ARITH_FOLDER) -DRRLWR_SECURITY_LEVEL=128 $(CFLAGS) $(SOURCES) $(UTILS_FOLDER)/drng.c test/cpucycles.c test/speed_print.c test/test_speed_KEM.c -o $@

test/test_speed_KEM256: $(SOURCES) $(HEADERS) test/cpucycles.h test/cpucycles.c test/speed_print.c test/speed_print.h $(UTILS_FOLDER)/drng.c $(UTILS_FOLDER)/drng.h test/test_speed_KEM.c
	$(CC) -I$(BASE_FOLDER) -I$(UTILS_FOLDER) -I$(ARITH_FOLDER) -DRRLWR_SECURITY_LEVEL=256 $(CFLAGS) $(SOURCES) $(UTILS_FOLDER)/drng.c test/cpucycles.c test/speed_print.c test/test_speed_KEM.c -o $@

test/test_speed_KEM512: $(SOURCES) $(HEADERS) test/cpucycles.h test/cpucycles.c test/speed_print.c test/speed_print.h $(UTILS_FOLDER)/drng.c $(UTILS_FOLDER)/drng.h test/test_speed_KEM.c
	$(CC) -I$(BASE_FOLDER) -I$(UTILS_FOLDER) -I$(ARITH_FOLDER) -DRRLWR_SECURITY_LEVEL=512 $(CFLAGS) $(SOURCES) $(UTILS_FOLDER)/drng.c test/cpucycles.c test/speed_print.c test/test_speed_KEM.c -o $@

KAT/KAT_KEM128: $(SOURCES) $(HEADERS) $(UTILS_FOLDER)/drng.c $(UTILS_FOLDER)/drng.h KAT/KAT_KEM.c
	$(CC) -I$(BASE_FOLDER) -I$(UTILS_FOLDER) -I$(ARITH_FOLDER) -DRRLWR_SECURITY_LEVEL=128 $(CFLAGS) $(SOURCES) $(UTILS_FOLDER)/drng.c KAT/KAT_KEM.c -o $@

KAT/KAT_KEM256: $(SOURCES) $(HEADERS) $(UTILS_FOLDER)/drng.c $(UTILS_FOLDER)/drng.h KAT/KAT_KEM.c
	$(CC) -I$(BASE_FOLDER) -I$(UTILS_FOLDER) -I$(ARITH_FOLDER) -DRRLWR_SECURITY_LEVEL=256 $(CFLAGS) $(SOURCES) $(UTILS_FOLDER)/drng.c KAT/KAT_KEM.c -o $@

KAT/KAT_KEM512: $(SOURCES) $(HEADERS) $(UTILS_FOLDER)/drng.c $(UTILS_FOLDER)/drng.h KAT/KAT_KEM.c
	$(CC) -I$(BASE_FOLDER) -I$(UTILS_FOLDER) -I$(ARITH_FOLDER) -DRRLWR_SECURITY_LEVEL=512 $(CFLAGS) $(SOURCES) $(UTILS_FOLDER)/drng.c KAT/KAT_KEM.c -o $@

clean:
	-$(RM) -f test/test_KEM128
	-$(RM) -f test/test_KEM256
	-$(RM) -f test/test_KEM512
	-$(RM) -f test/unit_tests_KEM128
	-$(RM) -f test/unit_tests_KEM256
	-$(RM) -f test/unit_tests_KEM512
	-$(RM) -f test/test_speed_KEM128
	-$(RM) -f test/test_speed_KEM256
	-$(RM) -f test/test_speed_KEM512
	-$(RM) -f KAT/KAT_KEM128
	-$(RM) -f KAT/KAT_KEM256
	-$(RM) -f KAT/KAT_KEM512
	-$(RM) -Rf output
