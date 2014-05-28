package io;

import it.unimi.dsi.fastutil.longs.Long2IntMap;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.dna.kmers.Kmer;
import ru.ifmo.genetics.dna.kmers.KmerIteratorFactory;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.executors.BlockingThreadPoolExecutor;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;

import java.io.*;
import java.util.Random;
import java.util.concurrent.CountDownLatch;

import org.apache.log4j.Logger;
import ru.ifmo.genetics.tools.io.LazyBinqReader;
import ru.ifmo.genetics.tools.io.LazyDnaReader;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;

/**
 * Vladimir Ulyantsev
 * Date: 11.02.14
 * Time: 17:32
 */
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
            BufferedReader reader = new BufferedReader(new FileReader(file));

            int cnt = 0;
            long len = 0;
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.charAt(0) != '>' && line.length() >= minSeqLen) {
                    addSequence(hm, line, k);
                    cnt++;
                    len += line.length();
                }
            }
            logger.debug(cnt + " sequences added, summary len = " + len + " from " + file.getPath());
            totalSeq += cnt;
            totalLen += len;
        }

        logger.debug("Total sequences count = " + totalSeq);
        logger.debug("Total sequences length = " + totalLen);
        logger.debug("k-mers HM size = " + hm.size());
    }

    private static void addSequence(ArrayLong2IntHashMap hm, String seq, int k) {
        ShortKmer kmer = new ShortKmer(seq.substring(0, k));
        hm.add(kmer.toLong(), 1);
        for (int pos = k; pos < seq.length(); pos++) {
            kmer.shiftRight(DnaTools.fromChar(seq.charAt(pos)));
            hm.add(kmer.toLong(), 1);
        }
    }

    public static void printKmers(ArrayLong2IntHashMap hm,
                                  String path,
                                  int threshold) throws IOException {
        DataOutputStream stream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(new File(path))));
        //DataOutputStream stream = new DataOutputStream(new FileOutputStream(new File(path)));

        for (Long2IntMap map : hm.hm) {
            for (Long2IntMap.Entry e : map.long2IntEntrySet()) {
                if (e.getIntValue() > threshold) {
                    stream.writeLong(e.getLongKey());
                    stream.writeInt(e.getIntValue());
                }
            }
        }

        stream.close();
    }

    public static ArrayLong2IntHashMap loadKmers(File[] files,
                                                 int freqThreshold,
                                                 int availableProcessors,
                                                 Logger logger) throws IOException {
        ArrayLong2IntHashMap hm =
                new ArrayLong2IntHashMap((int) (Math.log(availableProcessors) / Math.log(2)) + 4);
        addKmers(files, hm, freqThreshold, logger);
        return hm;
    }

    public static void addKmers(File[] files,
                                ArrayLong2IntHashMap hm,
                                int freqThreshold,
                                Logger logger) throws IOException {
        for (File file : files) {
            DataInputStream inputStream = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));

            long uniqueKmers = 0, uniqueKmersAdded = 0, totalKmers = 0, totalKmersAdded = 0;
            while (inputStream.available() > 2) {
                long kmerRepr = inputStream.readLong();
                int freq = inputStream.readInt();

                uniqueKmers++;
                totalKmers += freq;
                if (freq > freqThreshold) {
                    uniqueKmersAdded++;
                    totalKmersAdded += freq;
                    hm.add(kmerRepr, freq);
                }
            }

            inputStream.close();

            logger.debug(file + " : " + uniqueKmersAdded + " / " + uniqueKmers + " unique k-mers added");
            logger.debug(file + " : " + totalKmersAdded + " / " + totalKmers + " total k-mers added");
        }
        logger.debug("k-mers loaded");
    }

    public static ArrayLong2IntHashMap loadKmers(File[] files,
                                                 int freqThreshold,
                                                 int availableProcessors) throws InterruptedException {

        ArrayLong2IntHashMap hm =
                new ArrayLong2IntHashMap((int) (Math.log(availableProcessors) / Math.log(2)) + 4);
        addKmers(files, hm, freqThreshold, availableProcessors);
        return hm;
    }

    public static void addKmers(File[] files,
                                ArrayLong2IntHashMap hm,
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

    public static ArrayLong2IntHashMap loadBINQReads(File[] files,
                                                     int k,
                                                     int loadTaskSize,
                                                     KmerIteratorFactory<? extends Kmer> factory,
                                                     int availableProcessors,
                                                     Logger logger) throws IOException {
        ArrayLong2IntHashMap hm =
                new ArrayLong2IntHashMap((int) (Math.log(availableProcessors) / Math.log(2)) + 4);
        addBINQReads(files, hm, k, loadTaskSize, factory, availableProcessors, logger);
        return hm;
    }

    public static void addBINQReads(File[] files,
                                    ArrayLong2IntHashMap hm,
                                    int k,
                                    int loadTaskSize,
                                    KmerIteratorFactory<? extends Kmer> factory,
                                    int availableProcessors,
                                    Logger logger) throws IOException {
        LazyBinqReader reader = new LazyBinqReader(files);

        DnaQReadDispatcher dispatcher = new DnaQReadDispatcher(reader, loadTaskSize);
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
        logger.info(dispatcher.reads + " reads added");
        logger.info("k-mers loaded");
    }

    public static ArrayLong2IntHashMap loadReads(File[] files,
                                                 int k,
                                                 int loadTaskSize,
                                                 KmerIteratorFactory<? extends Kmer> factory,
                                                 int availableProcessors,
                                                 Logger logger) throws IOException, ExecutionFailedException {
        ArrayLong2IntHashMap hm =
                new ArrayLong2IntHashMap((int) (Math.log(availableProcessors) / Math.log(2)) + 4);
        addReads(files, hm, k, loadTaskSize, factory, availableProcessors, logger);
        return hm;
    }

    public static void addReads(File[] files,
                                ArrayLong2IntHashMap hm,
                                int k,
                                int loadTaskSize,
                                KmerIteratorFactory<? extends Kmer> factory,
                                int availableProcessors,
                                Logger logger) throws ExecutionFailedException {
        for (File file : files) {
            NamedSource<Dna> reader = LazyDnaReader.sourceFromFile(file);

            UniversalReadDispatcher dispatcher = new UniversalReadDispatcher(reader, loadTaskSize);
            UniversalLoadWorker[] workers = new UniversalLoadWorker[availableProcessors];
            CountDownLatch latch = new CountDownLatch(workers.length);

            for (int i = 0; i < workers.length; ++i) {
                workers[i] = new UniversalLoadWorker(dispatcher, latch, k, hm, factory);
                new Thread(workers[i]).start();
            }

            try {
                latch.await();
            } catch (InterruptedException e) {
                logger.warn("Main thread interrupted");
                for (UniversalLoadWorker worker : workers) {
                    worker.interrupt();
                }
            }
            logger.info(dispatcher.reads + " reads added from " + file);
        }
        logger.info("k-mers loaded");
    }

}
