CC ?= gcc
SECURITY_LEVEL ?= 128
ROOT := ..

CFLAGS ?= -O3
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -fomit-frame-pointer

INCLUDES := -I$(ROOT) -I$(ROOT)/arith -I$(ROOT)/utils -I$(ROOT)/test

COMMON_SOURCES := \
  $(ROOT)/arith/fprime.c \
  $(ROOT)/arith/poly.c \
  $(ROOT)/arith/ring.c \
  $(ROOT)/arith/packing.c \
  $(ROOT)/arith/uniform.c \
  $(ROOT)/pke.c \
  $(ROOT)/utils/fips202.c \
  $(ROOT)/utils/drng.c \
  $(ROOT)/test/cpucycles.c \
  $(ROOT)/test/speed_print.c

SOURCES := \
  $(COMMON_SOURCES) \
  matrixVectorMul_test.c

TARGET := matrixVectorMul_test$(SECURITY_LEVEL)
POLY_TARGET := polyMul_test$(SECURITY_LEVEL)

.PHONY: all run poly run-poly clean

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CC) $(INCLUDES) -DRRLWR_SECURITY_LEVEL=$(SECURITY_LEVEL) $(CFLAGS) $(SOURCES) -o $@

$(POLY_TARGET): $(COMMON_SOURCES) polyMul_test.c
	$(CC) $(INCLUDES) -DRRLWR_SECURITY_LEVEL=$(SECURITY_LEVEL) $(CFLAGS) $(COMMON_SOURCES) polyMul_test.c -o $@

run: $(TARGET)
	./$(TARGET)

poly: $(POLY_TARGET)

run-poly: $(POLY_TARGET)
	./$(POLY_TARGET)

clean:
	$(RM) matrixVectorMul_test128 matrixVectorMul_test256 matrixVectorMul_test512 \
	      polyMul_test128 polyMul_test256 polyMul_test512
