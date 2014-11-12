package io;

import it.unimi.dsi.fastutil.longs.Long2IntMap;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.dna.kmers.Kmer;
import ru.ifmo.genetics.dna.kmers.KmerIteratorFactory;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.dna.kmers.ShortKmerIteratorFactory;
import ru.ifmo.genetics.executors.BlockingThreadPoolExecutor;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.io.sources.Source;
import ru.ifmo.genetics.statistics.QuickQuantitativeStatistics;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;

import java.io.*;
import java.util.Iterator;
import java.util.Random;
import java.util.concurrent.CountDownLatch;

import org.apache.log4j.Logger;
import ru.ifmo.genetics.structures.map.BigLong2IntHashMap;
import ru.ifmo.genetics.structures.map.MutableLongIntEntry;
import ru.ifmo.genetics.tools.ec.DnaQReadDispatcher;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.Tool;

public class IOUtils {
    static final int BYTES_PER_KMER = 12;

    public static void addFASTASequences(File[] files,
                                         ArrayLong2IntHashMap hm,
                                         int k,
                                         int minSeqLen,
                                         Logger logger) throws IOException {
        int totalSeq = 0;
        long totalLen = 0;
        for (File file : files) {
            Source<Dna> source = ReadersUtils.readDnaLazy(file);
            Iterator<Dna> in = source.iterator();

            int cnt = 0;
            long len = 0;
            while (in.hasNext()) {
                Dna dna = in.next();
                if (dna.length() >= minSeqLen) {
                    addSequence(hm, dna, k);
                    cnt++;
                    len += dna.length();
                }
            }
            logger.debug(cnt + " sequences added, summary len = " + len + " from " + file.getName());
            totalSeq += cnt;
            totalLen += len;
        }

        logger.debug("Total sequences count = " + totalSeq);
        logger.debug("Total sequences length = " + totalLen);
        logger.debug("k-mers HM size = " + hm.size());
    }

    private static void addSequence(ArrayLong2IntHashMap hm, Dna dna, int k) {
        Iterator<ShortKmer> it = ShortKmer.kmersOf(dna, k).iterator();
        while (it.hasNext()) {
            hm.add(it.next().toLong(), 1);
        }
    }

    public static long printKmers(BigLong2IntHashMap hm, int threshold,
                                  File outFile, File stFile) throws IOException {
        DataOutputStream stream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile)));

        QuickQuantitativeStatistics<Integer> stats = new QuickQuantitativeStatistics<Integer>();
        long good = 0;

        Iterator<MutableLongIntEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongIntEntry entry = it.next();
            long key = entry.getKey();
            int value = entry.getValue();

            stats.add(value);

            if (value > threshold) {
                stream.writeLong(key);
                stream.writeInt(value);
                good++;
            }
        }

        stream.close();
        stats.printToFile(stFile, "# k-mer frequency\tnumber of such k-mers");
        return good;
    }

    public static BigLong2IntHashMap loadKmers(File[] files,
                                                 int freqThreshold,
                                                 int availableProcessors,
                                                 Logger logger) throws IOException {
        BigLong2IntHashMap hm = new BigLong2IntHashMap(
                (int) (Math.log(availableProcessors) / Math.log(2)) + 4, 8);
        addKmers(files, hm, freqThreshold, logger);
        return hm;
    }

    public static void addKmers(File[] files,
                                BigLong2IntHashMap hm,
                                int freqThreshold,
                                Logger logger) throws IOException {
        for (File file : files) {
            Tool.info(logger, "Loading k-mer from " + file + "...");

            FileInputStream fis = new FileInputStream(file);
            DataInputStream is = new DataInputStream(new BufferedInputStream(fis));

            long uniqueKmers = 0, uniqueKmersAdded = 0, totalKmers = 0, totalKmersAdded = 0;
            long c = fis.getChannel().size() / 12;
            uniqueKmers = c;
            for (; c > 0; c--) {
                long kmerRepr = is.readLong();
                int freq = is.readInt();

//                uniqueKmers++;
                totalKmers += freq;
                if (freq > freqThreshold) {
                    uniqueKmersAdded++;
                    totalKmersAdded += freq;
                    hm.addAndBound(kmerRepr, freq);
                }
            }

            if (is.available() > 2) {
                throw new RuntimeException("Size mismatch. Possibly wrong file format/file is corrupted.");
            }
            is.close();

            Tool.debug(logger, file + " : " + NumUtils.groupDigits(uniqueKmersAdded) + " / "
                    + NumUtils.groupDigits(uniqueKmers) + " unique k-mers added/all");
            Tool.debug(logger, file + " : " + NumUtils.groupDigits(totalKmersAdded) + " / "
                    + NumUtils.groupDigits(totalKmers) + " total k-mers added/all");
            Tool.info(logger, NumUtils.groupDigits(uniqueKmersAdded) + " k-mers loaded from " + file);
        }
    }

    public static BigLong2IntHashMap loadKmers(File[] files,
                                                 int freqThreshold,
                                                 int availableProcessors) throws InterruptedException {

        BigLong2IntHashMap hm = new BigLong2IntHashMap(
                (int) (Math.log(availableProcessors) / Math.log(2)) + 4, 8);
        addKmers(files, hm, freqThreshold, availableProcessors);
        return hm;
    }

    public static void addKmers(File[] files,
                                BigLong2IntHashMap hm,
                                int freqThreshold,
                                int availableProcessors) throws InterruptedException {

        BlockingThreadPoolExecutor executor = new BlockingThreadPoolExecutor(availableProcessors);

        for (File file : files) {
            long bytesInFile = file.length();
            long kmersCount = bytesInFile / BYTES_PER_KMER;
            long kmersPerThread = kmersCount / availableProcessors + 1;

            long kmersSum = 0;
            for (int i = 0; i < availableProcessors; i++) {
                long kmersToAdd = Math.min(kmersPerThread, kmersCount - kmersSum);
                executor.blockingExecute(new KmersToHMAdditionTask(hm, file, kmersSum, kmersToAdd, freqThreshold));
                kmersSum += kmersToAdd;
            }
        }

        executor.shutdownAndAwaitTermination();
    }

    public static BigLong2IntHashMap loadBINQReads(File[] files,
                                                     int k,
                                                     int loadTaskSize,
                                                     KmerIteratorFactory<? extends Kmer> factory,
                                                     int availableProcessors,
                                                     Logger logger) throws IOException {
        BigLong2IntHashMap hm = new BigLong2IntHashMap(
                (int) (Math.log(availableProcessors) / Math.log(2)) + 4, 8);
        addBINQReads(files, hm, k, loadTaskSize, factory, availableProcessors, logger);
        return hm;
    }

    public static void addBINQReads(File[] files,
                                    BigLong2IntHashMap hm,
                                    int k,
                                    int loadTaskSize,
                                    KmerIteratorFactory<? extends Kmer> factory,
                                    int availableProcessors,
                                    Logger logger) throws IOException {
        Source<DnaQ> source = ReadersUtils.readDnaQLazy(files);
        Iterator<DnaQ> it = source.iterator();


        DnaQReadDispatcher dispatcher = new DnaQReadDispatcher(source, loadTaskSize, null);
        KmerLoadWorker[] workers = new KmerLoadWorker[availableProcessors];
        CountDownLatch latch = new CountDownLatch(workers.length);

        for (int i = 0; i < workers.length; ++i) {
            workers[i] = new KmerLoadWorker(dispatcher, latch, new Random(42),
                    k, hm, factory);
            new Thread(workers[i]).start();
        }

        try {
            latch.await();
        } catch (InterruptedException e) {
            logger.warn("Main thread interrupted");
            for (KmerLoadWorker worker : workers) {
                worker.interrupt();
            }
        }
        logger.info(NumUtils.groupDigits(dispatcher.getReads()) + " reads added");
//        logger.info("k-mers loaded");
    }

    public static BigLong2IntHashMap loadReads(File[] files,
                                                 int k,
                                                 int loadTaskSize,
                                                 KmerIteratorFactory<? extends Kmer> factory,
                                                 int availableProcessors,
                                                 Logger logger) throws IOException, ExecutionFailedException {
        BigLong2IntHashMap hm = new BigLong2IntHashMap(
                (int) (Math.log(availableProcessors) / Math.log(2)) + 4, 8);

        addReads(files, hm, k, loadTaskSize, factory, availableProcessors, logger);
        return hm;
    }

    public static void addReads(File[] files,
                                BigLong2IntHashMap hm,
                                int k,
                                int loadTaskSize,
                                KmerIteratorFactory<? extends Kmer> factory,
                                int availableProcessors,
                                Logger logger) throws ExecutionFailedException {
        for (File file : files) {
            logger.info("Processing file " + file + "...");

            NamedSource<Dna> reader = ReadersUtils.readDnaLazy(file);

            UniversalReadDispatcher dispatcher = new UniversalReadDispatcher(reader, loadTaskSize, hm);
            UniversalLoadWorker[] workers = new UniversalLoadWorker[availableProcessors];
            CountDownLatch latch = new CountDownLatch(workers.length);

            logger.debug("Starting workers...");
            for (int i = 0; i < workers.length; ++i) {
                workers[i] = new UniversalLoadWorker(dispatcher, latch, k, hm, factory);
                new Thread(workers[i]).start();
            }


            try {
                logger.debug("Waiting workers...");
                latch.await();
            } catch (InterruptedException e) {
                logger.warn("Main thread interrupted");
                for (UniversalLoadWorker worker : workers) {
                    worker.interrupt();
                }
            }
            logger.info(NumUtils.groupDigits(dispatcher.reads) + " reads added from " + file);
        }
//        logger.info("k-mers loaded");
    }

}
