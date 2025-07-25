CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -g
INCLUDES = -I../../vendor/ITU-R-HF/P533/Src/P533 -I../../vendor/cjson
LIBS = -lm -L../../vendor/ITU-R-HF/ITURHFProp/Linux -lp533 -lp372
LDFLAGS = -Wl,-rpath,../../vendor/ITU-R-HF/ITURHFProp/Linux

SOURCES = c_reference_test.c
JSON_SOURCES = c_reference_test.c ../../vendor/cjson/cJSON.c
TARGET = c_reference_test
JSON_TARGET = c_reference_test

$(TARGET): $(SOURCES)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $(TARGET) $(SOURCES) $(LIBS)

$(JSON_TARGET): $(JSON_SOURCES)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $(JSON_TARGET) $(JSON_SOURCES) $(LIBS)

all: $(TARGET) $(JSON_TARGET)

clean:
	rm -f $(TARGET) $(JSON_TARGET) *.json *.o

test-json: $(JSON_TARGET)
	./$(JSON_TARGET)

.PHONY: clean test test-json all
