package io;

import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import ru.ifmo.genetics.structures.map.BigLong2IntHashMap;

import java.io.*;

/**
 * Created by ulyantsev on 05.05.14.
 *
 */
public class KmersToHMAdditionTask implements Runnable {

    final int BYTES_PER_KMER = 12;

    BigLong2IntHashMap hm;
    File kmersFile;
    long kmersToSkip;
    long kmersToRead;
    int maxBadFrequency;

    public KmersToHMAdditionTask(BigLong2IntHashMap hm,
                                 File kmersFile,
                                 long kmersToSkip,
                                 long kmersToRead,
                                 int maxBadFrequency) {
        this.hm = hm;
        this.kmersFile = kmersFile;
        this.kmersToSkip = kmersToSkip;
        this.kmersToRead = kmersToRead;
        this.maxBadFrequency = maxBadFrequency;
    }

    @Override
    public void run() {
        try {
            DataInputStream inputStream = new DataInputStream(new BufferedInputStream(new FileInputStream(kmersFile)));
            long reallyReaden = inputStream.skip(BYTES_PER_KMER * kmersToSkip);
            assert reallyReaden == BYTES_PER_KMER * kmersToSkip;

            for (long i = 0; i < kmersToRead; i++) {
                long kmerRepr = inputStream.readLong();
                int frequency = inputStream.readInt();
                if (frequency > maxBadFrequency) {
                    hm.addAndBound(kmerRepr, frequency);
                }
            }

            inputStream.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
