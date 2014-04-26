package io;

import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.io.IOUtils;
import ru.ifmo.genetics.tools.io.LazyBinqReader;
import ru.ifmo.genetics.utils.tool.Progress;

import java.io.EOFException;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

public class DnaQReadDispatcher {
    LazyBinqReader reader;
    int workRangeSize;
    final OutputStream out;

    long reads = 0;
    long totalReadsProcessed;

    long[] ar;
    int r = 0;
    int w = 0;

    public DnaQReadDispatcher(LazyBinqReader reader, int workRangeSize) {
        this.reader = reader;
        this.workRangeSize = workRangeSize;
        out = null;
    }

    public DnaQReadDispatcher(LazyBinqReader reader, OutputStream out, int workRangeSize, long totalReadsProcessed,
                              int workersNumber) {
        this.reader = reader;
        this.workRangeSize = workRangeSize;
        this.out = out;
        this.totalReadsProcessed = totalReadsProcessed;


        ar = new long[workersNumber];
    }

    private List<DnaQ> readRange(LazyBinqReader reader, int workRangeSize) throws IOException {
        List<DnaQ> list = new ArrayList<DnaQ>();
        while (list.size() < workRangeSize) {
            try {
                list.add(reader.readDnaq());
                ++reads;
            } catch (EOFException e) {
                break;
            }
        }
        return list;
    }


    public List<DnaQ> getWorkRange(long id) {
        List<DnaQ> list;
        try {
            synchronized (reader) {
                list = readRange(reader, workRangeSize);
                ar[r] = id;
                r = (r + 1) % ar.length;
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return list.isEmpty() ? null : list;
    }

    public List<DnaQ> getWorkRange() {
        List<DnaQ> list;
        try {
            synchronized (reader) {
                list = readRange(reader, workRangeSize);
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
                IOUtils.putByteArray(dnaq.toByteArray(), out);
                ++totalReadsProcessed;
            }
            w = (w + 1) % ar.length;
            out.notifyAll();
        }
    }

    public long getReads(){
        return reads;
    }

}
