package io;

import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.statistics.QuickQuantitativeStatistics;
import ru.ifmo.genetics.structures.map.BigLong2LongHashMap;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Tool;

import java.io.*;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.CountDownLatch;

public class IOUtils {

    static final int READS_WORK_RANGE_SIZE = 1 << 15;   // 32 K reads
    static final int KMERS_WORK_RANGE_SIZE = 16777220;   // ~16 Mb of data



    public static String withP(long cur, long all) {
        return NumUtils.groupDigits(cur) + " (" + String.format("%.1f", cur * 100.0 / all) + "%)";
    }
    public static String withP(long cur, long all, String firstAddition, String secondAddition) {
        return NumUtils.groupDigits(cur) + " " + firstAddition + " " +
                "(" + String.format("%.1f", cur * 100.0 / all) + "% " + secondAddition + ")";
    }




    public static long printKmers(BigLong2ShortHashMap hm, int threshold,
                                  File outFile, File stFile) throws IOException {
        DataOutputStream stream = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(outFile), 1 << 24));   // 16 Mb buffer

        QuickQuantitativeStatistics<Short> stats = new QuickQuantitativeStatistics<Short>();
        long good = 0;

        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            long key = entry.getKey();
            short value = entry.getValue();

            stats.add(value);

            if (value > threshold) {
                stream.writeLong(key);
                stream.writeShort(value);
                good++;
            }
        }

        stream.close();
        stats.printToFile(stFile, "# k-mer frequency\tnumber of such k-mers");
        return good;
    }

    public static long filterAndPrintKmers(BigLong2ShortHashMap hm, BigLong2ShortHashMap filter_hm,
                                           int threshold, int filter_threshold, File out) throws IOException {
        DataOutputStream stream = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(out), 1 << 24));   // 16 Mb buffer

        long good = 0;

        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            long key = entry.getKey();
            short value = entry.getValue();

            if (value > threshold && filter_hm.getWithZero(key) > filter_threshold) {
                stream.writeLong(key);
                stream.writeShort(value);
                good++;
            }
        }

        stream.close();
        return good;
    }

    public static long MultipleFiltersAndPrintKmers(BigLong2ShortHashMap hm,
                                                    BigLong2ShortHashMap cd_filter_hm,
                                                    BigLong2ShortHashMap uc_filter_hm,
                                                    BigLong2ShortHashMap nonibd_filter_hm,
                                                    int threshold, File out, File stFile) throws IOException {
        DataOutputStream stream = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(out), 1 << 24));   // 16 Mb buffer

        QuickQuantitativeStatistics<Triple> stats = new QuickQuantitativeStatistics<Triple>();
        long good = 0;

        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            long key = entry.getKey();
            short value = entry.getValue();

            if (value > threshold) {
                short cd = cd_filter_hm.getWithZero(key);
                short uc = uc_filter_hm.getWithZero(key);
                short nonibd = nonibd_filter_hm.getWithZero(key);

                stats.add(new Triple(cd, uc, nonibd));

                if (cd > 0 || uc > 0 || nonibd > 0) {
                    stream.writeLong(key);
                    stream.writeShort(value);
                    good++;
                }
            }
        }

        stream.close();
        stats.printToFile(stFile, "# cd k-mer samples\t" +
                "uc k-mer samples\tnonIBD k-mer samples\tnumber of such k-mers");
        return good;
    }

    private static class Triple implements Comparable<Triple> {
        private final short cd;
        private final short uc;
        private final short nonibd;

        private Triple(short cd, short uc, short nonibd) {
            this.cd = cd;
            this.uc = uc;
            this.nonibd = nonibd;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Triple)) return false;
            Triple triple = (Triple) o;
            return cd == triple.cd && uc == triple.uc && nonibd == triple.nonibd;
        }

        @Override
        public int hashCode() {
            int result = cd;
            result = 31 * result + uc;
            result = 31 * result + nonibd;
            return result;
        }

        @Override
        public int compareTo(Triple o) {
            if (this.equals(o)) return 0;

            if (this.cd < o.cd) return -1;
            else {
                if (this.cd > o.cd) return 1;
                else {
                    if (this.uc < o.uc) return -1;
                    else {
                        if (this.uc > o.uc) return 1;
                        else {
                            return Short.compare(this.nonibd, o.nonibd);
                        }
                    }
                }
            }
        }

        @Override
        public String toString() {
            return cd + "\t" + uc + "\t" + nonibd;
        }
    }



    public static void tryToAppendDescription(File[] outputFilesDesc, File f, String msg) {
        if (outputFilesDesc != null) {
            for (File wf : outputFilesDesc) {
                try {
                    PrintWriter out = new PrintWriter(new FileWriter(wf, true));
                    out.println();
                    out.println(f);
                    out.println("   " + msg);
                    out.close();
                } catch (IOException e) {
                    // does not matter
                }
            }
        }
    }



    // ---------------------------- for loading kmers ----------------------------------

    static class Kmers2HMWorker extends KmersLoadWorker {
        Kmers2HMWorker(BigLong2ShortHashMap hm, int freqThreshold) {
            this.hm = hm;
            this.freqThreshold = freqThreshold;
        }

        final BigLong2ShortHashMap hm;
        final int freqThreshold;
        long kmers = 0, kmersAdded = 0;
        long freqSum = 0, freqSumAdded = 0;

        @Override
        public void processKmer(long kmer, short freq) {
            kmers++;
            freqSum += freq;
            if (freq > freqThreshold) {
                hm.addAndBound(kmer, freq);
                kmersAdded++;
                freqSumAdded += freq;
            }
        }
    }

    public static BigLong2ShortHashMap loadKmers(File[] files, int freqThreshold, int availableProcessors, Logger logger)
            throws ExecutionFailedException {

        BigLong2ShortHashMap hm = new BigLong2ShortHashMap(
                (int) (Math.log(availableProcessors) / Math.log(2)) + 4, 12);

        Kmers2HMWorker[] workers = new Kmers2HMWorker[availableProcessors];
        for (int i = 0; i < workers.length; ++i) {
            workers[i] = new Kmers2HMWorker(hm, freqThreshold);
        }

        run(files, workers, hm, logger);

        // calculating statistics...
        long kmers = 0, kmersAdded = 0;
        long freqSum = 0, freqSumAdded = 0;
        for (Kmers2HMWorker worker : workers) {
            kmers += worker.kmers;
            kmersAdded += worker.kmersAdded;
            freqSum += worker.freqSum;
            freqSumAdded += worker.freqSumAdded;
        }

        Tool.debug(logger,
                "Added/All kmers count = " + NumUtils.groupDigits(kmersAdded) + "/" + NumUtils.groupDigits(kmers)
                        + " (" + String.format("%.1f", kmersAdded * 100.0 / kmers) + "%)");
        Tool.debug(logger,
                "Added/All kmers frequency sum = " + NumUtils.groupDigits(freqSumAdded) + "/" + NumUtils.groupDigits(freqSum)
                        + " (" + String.format("%.1f", freqSumAdded * 100.0 / freqSum) + "%)");
        Tool.debug(logger, "k-mers HM size = " + NumUtils.groupDigits(hm.size()));

        return hm;
    }


    static class KmersPresenceWorker extends KmersLoadWorker {
        KmersPresenceWorker(BigLong2LongHashMap hm) {
            this.hm = hm;
        }
        final BigLong2LongHashMap hm;
        @Override
        public void processKmer(long kmer, short freq) {
            if (hm.contains(kmer)) {
                hm.addAndBound(kmer, freq);
            }
        }
    }

    public static void calculatePresenceForKmers(File[] files, BigLong2LongHashMap hm, int availableProcessors, Logger logger)
            throws ExecutionFailedException {
        BytesWorker[] workers = new BytesWorker[availableProcessors];
        for (int i = 0; i < workers.length; ++i) {
            workers[i] = new KmersPresenceWorker(hm);
        }
        run(files, workers, null, logger);
    }


    public static void run(File[] files, BytesWorker[] workers, BigLong2ShortHashMap hmForMonitoring, Logger logger)
            throws ExecutionFailedException {
        try {
            for (File file : files) {
                Tool.info(logger, "Loading file " + file.getName() + "...");

                InputStream is = new FileInputStream(file);
                BytesDispatcher dispatcher = new BytesDispatcher(is, KMERS_WORK_RANGE_SIZE, hmForMonitoring);
                CountDownLatch latch = new CountDownLatch(workers.length);

                for (int i = 0; i < workers.length; ++i) {
                    workers[i].setDispatcher(dispatcher);
                    workers[i].setLatch(latch);
                    new Thread(workers[i]).start();
                }

                try {
                    latch.await();
                } catch (InterruptedException e) {
                    Tool.warn(logger, "Main thread interrupted");
                    for (BytesWorker worker : workers) {
                        worker.interrupt();
                    }
                    throw new ExecutionFailedException("Thread was interrupted", e);
                }
                Tool.debug(logger, NumUtils.memoryAsString(dispatcher.bytesRead) + " of data processed");
            }
        } catch (IOException e) {
            throw new ExecutionFailedException("Can't load k-mers file", e);
        }
    }


    // ---------------------------- for loading reads ----------------------------------

    static class ReadsLoadWorker extends ReadsWorker {
        ReadsLoadWorker(BigLong2ShortHashMap hm, int k, int minDnaLen) {
            this.hm = hm;
            this.k = k;
            this.minDnaLen = minDnaLen;
        }

        final BigLong2ShortHashMap hm;
        final int k;
        final int minDnaLen;
        int totalSeq = 0, goodSeq = 0;
        long totalLen = 0, goodLen = 0;

        @Override
        public void process(List<Dna> reads) {
            for (Dna dna : reads) {
                totalSeq++;
                totalLen += dna.length();

                if (dna.length() >= minDnaLen) {
                    for (ShortKmer kmer : ShortKmer.kmersOf(dna, k)) {
                        hm.addAndBound(kmer.toLong(), (short) 1);
                    }
                    goodSeq++;
                    goodLen += dna.length();
                }
            }
        }
    }

    public static BigLong2ShortHashMap loadReads(File[] files, int k, int minSeqLen,
                                                 int availableProcessors, Logger logger)
            throws ExecutionFailedException, IOException {
        BigLong2ShortHashMap hm = new BigLong2ShortHashMap(
                (int) (Math.log(availableProcessors) / Math.log(2)) + 4, 12, true);

        ReadsLoadWorker[] workers = new ReadsLoadWorker[availableProcessors];
        for (int i = 0; i < workers.length; ++i) {
            workers[i] = new ReadsLoadWorker(hm, k, minSeqLen);
        }

        run(files, workers, hm, logger);

        // calculating statistics...
        int totalSeq = 0, goodSeq = 0;
        long totalLen = 0, goodLen = 0;
        for (ReadsLoadWorker worker : workers) {
            totalSeq += worker.totalSeq;
            goodSeq += worker.goodSeq;
            totalLen += worker.totalLen;
            goodLen += worker.goodLen;
        }
        Tool.debug(logger,
                "Good/Total sequences count = " + NumUtils.groupDigits(goodSeq) + "/" + NumUtils.groupDigits(totalSeq)
                + " (" + String.format("%.1f", goodSeq * 100.0 / totalSeq) + "%)");
        Tool.debug(logger,
                "Good/Total sequences length = " + NumUtils.groupDigits(goodLen) + "/" + NumUtils.groupDigits(totalLen)
                        + " (" + String.format("%.1f", goodLen * 100.0 / totalLen) + "%)");
        Tool.debug(logger, "k-mers HM size = " + NumUtils.groupDigits(hm.size()));

        return hm;
    }


    static class ReadsPresenceWorker extends ReadsWorker {
        ReadsPresenceWorker(BigLong2LongHashMap hm, int k) {
            this.hm = hm;
            this.k = k;
        }

        final BigLong2LongHashMap hm;
        final int k;

        @Override
        public void process(List<Dna> reads) {
            for (Dna dna : reads) {
                for (ShortKmer kmer : ShortKmer.kmersOf(dna, k)) {
                    if (hm.contains(kmer.toLong())) {
                        hm.addAndBound(kmer.toLong(), 1);
                    }
                }
            }
        }
    }

    public static void calculatePresenceForReads(File[] files, int k, BigLong2LongHashMap hm, int availableProcessors, Logger logger)
            throws ExecutionFailedException, IOException {
        ReadsWorker[] workers = new ReadsWorker[availableProcessors];
        for (int i = 0; i < workers.length; ++i) {
            workers[i] = new ReadsPresenceWorker(hm, k);
        }
        run(files, workers, null, logger);
    }



    public static void run(File[] files, ReadsWorker[] workers, BigLong2ShortHashMap hmForMonitoring, Logger logger)
            throws ExecutionFailedException, IOException {
        for (File file : files) {
            Tool.info(logger, "Loading file " + file.getName() + "...");

            NamedSource<Dna> reader = ReadersUtils.readDnaLazy(file);

            ReadsDispatcher dispatcher = new ReadsDispatcher(reader, READS_WORK_RANGE_SIZE, hmForMonitoring);
            CountDownLatch latch = new CountDownLatch(workers.length);

            for (int i = 0; i < workers.length; ++i) {
                workers[i].setDispatcher(dispatcher);
                workers[i].setLatch(latch);
                new Thread(workers[i]).start();
            }

            try {
                latch.await();
            } catch (InterruptedException e) {
                Tool.warn(logger, "Main thread interrupted");
                for (ReadsWorker worker : workers) {
                    worker.interrupt();
                }
                throw new ExecutionFailedException("Thread was interrupted", e);
            }
            Tool.info(logger, NumUtils.groupDigits(dispatcher.reads) + " reads added");
        }
    }

}
