package io;

import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.dna.kmers.Kmer;
import ru.ifmo.genetics.dna.kmers.KmerIteratorFactory;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;

import java.util.List;
import java.util.Random;
import java.util.concurrent.CountDownLatch;

public class UniversalLoadWorker implements Runnable {

    private UniversalReadDispatcher dispatcher;
    CountDownLatch latch;
    int LEN;
    long step;
    ArrayLong2IntHashMap hm;
    KmerIteratorFactory<? extends Kmer> factory;

    boolean interrupted = false;

    Logger logger;

    public UniversalLoadWorker(UniversalReadDispatcher dispatcher, CountDownLatch latch,
                               int LEN, ArrayLong2IntHashMap hm,
                               KmerIteratorFactory<? extends Kmer> factory) {

        this.dispatcher = dispatcher;
        this.latch = latch;
        this.LEN = LEN;
        this.hm = hm;

        this.factory = factory;
    }



    void add(Iterable<Dna> dnas) {
        for (Dna dna : dnas) {
            for (ShortKmer kmer : ShortKmer.kmersOf(dna, LEN)) {
                hm.add(kmer.toLong(), 1);
            }
        }
    }

    public void interrupt() {
        interrupted = true;
    }

    public void run() {
        logger = Logger.getLogger("worker-" + Thread.currentThread().getId());
        while (!interrupted) {
            List<Dna> list = dispatcher.getWorkRange();
            if (list == null) {
                break;
            }

            add(list);
        }
        latch.countDown();
    }
}
