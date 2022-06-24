package io;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;

public abstract class LongKmersLoadWorker extends BytesWorker {

    final static int KMER_RECORD_SIZE = 16;

    public abstract void processKmer(long kmer, long freq);

    @Override
    public void process(byte[] range, int len) {
        if (len % KMER_RECORD_SIZE != 0) {
            throw new RuntimeException("BAD division by work range");
        }
        try {
            DataInputStream is = new DataInputStream(new ByteArrayInputStream(range, 0, len));
            int c = len / KMER_RECORD_SIZE;
            for (int i = 0; i < c; i++) {
                long kmer = is.readLong();
                long freq = is.readLong();
                processKmer(kmer, freq);
            }
            if (is.available() > 0) {
                throw new RuntimeException("Bad KMER_RECORD_SIZE!");
            }
        } catch (IOException e) {
            throw new RuntimeException("Can't load kmers from file", e);
        }
    }
}
