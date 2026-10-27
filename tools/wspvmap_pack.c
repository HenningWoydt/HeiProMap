/* Packs a raw little-endian u32 vertex -> PE array into a WSPVMAP1 file: the 72-byte header of a base file with the
 * payload hash (FNV-1a over the bytes of the payload, bytes 64-71) recomputed, followed by the payload.
 * Same as write_wspvmap1 in wse-spmv/run_all_mappers.py, which is too slow in Python for 10M entries.
 *
 * usage: wspvmap_pack BASE.wspvmap RAW.u32 OUT.wspvmap1
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define HEADER_BYTES 72

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s BASE RAW OUT\n", argv[0]);
        return 2;
    }
    unsigned char header[HEADER_BYTES];
    FILE *base = fopen(argv[1], "rb");
    if (!base || fread(header, 1, HEADER_BYTES, base) != HEADER_BYTES) {
        fprintf(stderr, "cannot read the header of %s\n", argv[1]);
        return 1;
    }
    fclose(base);

    FILE *raw = fopen(argv[2], "rb");
    if (!raw) {
        fprintf(stderr, "cannot open %s\n", argv[2]);
        return 1;
    }
    fseek(raw, 0, SEEK_END);
    long size = ftell(raw);
    fseek(raw, 0, SEEK_SET);
    unsigned char *payload = malloc((size_t) size);
    if (!payload || fread(payload, 1, (size_t) size, raw) != (size_t) size) {
        fprintf(stderr, "cannot read %s\n", argv[2]);
        return 1;
    }
    fclose(raw);

    uint64_t h = 14695981039346656037ULL;
    for (long i = 0; i < size; ++i) {
        h ^= payload[i];
        h *= 1099511628211ULL;
    }
    for (int i = 0; i < 8; ++i) { header[64 + i] = (unsigned char) (h >> (8 * i)); }

    FILE *out = fopen(argv[3], "wb");
    if (!out || fwrite(header, 1, HEADER_BYTES, out) != HEADER_BYTES
        || fwrite(payload, 1, (size_t) size, out) != (size_t) size) {
        fprintf(stderr, "cannot write %s\n", argv[3]);
        return 1;
    }
    fclose(out);
    free(payload);
    return 0;
}
