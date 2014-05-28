package io;

import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.utils.iterators.ProgressableIterator;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;

public class UniversalReadDispatcher {
    NamedSource<Dna> reader;
    ProgressableIterator<Dna> iterator;
    int workRangeSize;
    final OutputStream out;

    long reads = 0;
    long totalReadsProcessed;

    long[] ar;
    int r = 0;
    int w = 0;

    public UniversalReadDispatcher(NamedSource<Dna> reader, int workRangeSize) {
        this.reader = reader;
        this.iterator = reader.iterator();
        this.workRangeSize = workRangeSize;
        out = null;
    }

    public UniversalReadDispatcher(NamedSource<Dna> reader,
                                   OutputStream out,
                                   int workRangeSize,
                                   long totalReadsProcessed,
                                   int workersNumber) {
        this.reader = reader;
        this.iterator = reader.iterator();

        this.workRangeSize = workRangeSize;
        this.out = out;
        this.totalReadsProcessed = totalReadsProcessed;

        ar = new long[workersNumber];
    }

    private List<Dna> readRange(int workRangeSize) throws IOException {
        List<Dna> list = new ArrayList<Dna>();
        while (this.iterator != null && list.size() < workRangeSize) {
            try {
                list.add(iterator.next());
                ++reads;
            } catch (NoSuchElementException e) {
                this.iterator = null;
                break;
            } catch (IllegalArgumentException ignored) {

            }
//            //This part is needed due to buggy iterator.hasNext()
//            if (!this.iterator.hasNext()) {
//                this.iterator = null;
//                break;
//            }
//            list.add(iterator.next());
        }
        return list;
    }

    public List<Dna> getWorkRange(long id) {
        List<Dna> list;
        try {
            synchronized (reader) {
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
            synchronized (reader) {
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
