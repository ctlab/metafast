package io;

import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.io.sources.Source;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import ru.ifmo.genetics.structures.map.BigLong2IntHashMap;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.TextUtils;
import ru.ifmo.genetics.utils.iterators.ProgressableIterator;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;

public class UniversalReadDispatcher {
    protected final Logger logger = Logger.getLogger("read-dispatcher");

    ProgressableIterator<Dna> iterator;
    int workRangeSize;
    final OutputStream out;

    long reads = 0;
    long totalReadsProcessed;

    long[] ar;
    int r = 0;
    int w = 0;
    final BigLong2IntHashMap hm; // tmp

    public UniversalReadDispatcher(Source<Dna> reader, int workRangeSize, BigLong2IntHashMap hm) {
        this.iterator = reader.iterator();
        this.workRangeSize = workRangeSize;
        out = null;
        logger.debug("Using " + workRangeSize + " reads as workRangeSize");
        this.hm = hm;
    }

    public UniversalReadDispatcher(Source<Dna> reader,
                                   OutputStream out,
                                   int workRangeSize,
                                   long totalReadsProcessed,
                                   int workersNumber) {
        this.iterator = reader.iterator();

        this.workRangeSize = workRangeSize;
        this.out = out;
        this.totalReadsProcessed = totalReadsProcessed;

        ar = new long[workersNumber];
        hm = null;
    }

    private List<Dna> readRange(int workRangeSize) throws IOException {
        List<Dna> list = new ArrayList<Dna>(workRangeSize);
        while ((list.size() < workRangeSize) && iterator.hasNext()) {
            list.add(iterator.next());
            ++reads;

            if (reads % 1000000 == 0) {
                logger.debug("Processed " + NumUtils.groupDigits(reads) + " reads:");
                logger.debug("Total hm size = " + NumUtils.groupDigits(hm.size()) + ", " +
                        "size in hm.hm = {" + NumUtils.groupDigits(hm.maps[0].size()) + ", "
                        + NumUtils.groupDigits(hm.maps[1].size()) + ", "
                        + NumUtils.groupDigits(hm.maps[2].size()) + ", "
                        + NumUtils.groupDigits(hm.maps[3].size()) + ", ...}");
                logger.debug("Available memory (without running GC) = " + Misc.availableMemoryWithoutRunningGCAsString());
            }
        }
        return list;
    }

    public List<Dna> getWorkRange(long id) {
        List<Dna> list;
        try {
            synchronized (iterator) {
                list = readRange(workRangeSize);
                ar[r] = id;
                r = (r + 1) % ar.length;
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return list.isEmpty() ? null : list;
    }

    public List<Dna> getWorkRange() {
        List<Dna> list;
        try {
            synchronized (iterator) {
                list = readRange(workRangeSize);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return list.isEmpty() ? null : list;
    }

    public void writeDnaQs(List<DnaQ> list, long id) throws IOException, InterruptedException {
        synchronized (out) {
            while (ar[w] != id) {
                out.wait();
            }
            for (DnaQ dnaq : list) {
                ru.ifmo.genetics.io.IOUtils.putByteArray(dnaq.toByteArray(), out);
                ++totalReadsProcessed;
            }
            w = (w + 1) % ar.length;
            out.notifyAll();
        }
    }

    public long getReads() {
        return reads;
    }

}
