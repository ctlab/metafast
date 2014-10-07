package io;

import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.dna.kmers.*;
import ru.ifmo.genetics.dna.kmers.KmerIteratorFactory;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import ru.ifmo.genetics.tools.ec.DnaQReadDispatcher;

import java.util.List;
import java.util.Random;
import java.util.concurrent.CountDownLatch;

public class KmerLoadWorker implements Runnable {

    private DnaQReadDispatcher dispatcher;
    CountDownLatch latch;
    int LEN;
    long step;
    ArrayLong2IntHashMap hm;
    KmerIteratorFactory<? extends Kmer> factory;

    boolean interrupted = false;

    Random random;

    Logger logger;

    public KmerLoadWorker(DnaQReadDispatcher dispatcher, CountDownLatch latch, Random random,
                          int LEN, ArrayLong2IntHashMap hm,
                          KmerIteratorFactory<? extends Kmer> factory) {

        this.dispatcher = dispatcher;
        this.latch = latch;
        this.LEN = LEN;
        this.hm = hm;
        this.random = random;


        this.factory = factory;
    }

    public void add(Kmer kmer) {
        long key = kmer.toLong();
        hm.add(key, 1);
    }

    void add(DnaQ dnaq) {
        ShortKmer kmer = new ShortKmer(0, LEN);
        for (int pos = 0; pos < dnaq.length(); pos++) {
            kmer.shiftRight(dnaq.nucAt(pos));
            if (pos >= LEN - 1) {
                add(kmer);
            }
        }
    }

    void add(Iterable<DnaQ> dnaqs) {
        for (DnaQ dnaq : dnaqs) {
            add(dnaq);
        }
    }

    public void interrupt() {
        interrupted = true;
    }

    public void run() {
        logger = Logger.getLogger("worker-" + Thread.currentThread().getId());
        while (!interrupted) {
            List<DnaQ> list = dispatcher.getWorkRange();
            if (list == null) {
                break;
            }

            add(list);
        }
        latch.countDown();
    }

}
