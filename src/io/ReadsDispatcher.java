package io;

import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.DnaQ;
import ru.ifmo.genetics.io.sources.Source;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.iterators.ProgressableIterator;
import ru.ifmo.genetics.utils.tool.Tool;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class ReadsDispatcher {
    final Logger logger = Logger.getLogger("reads-dispatcher");

    final ProgressableIterator<Dna> iterator;
    public final int workRangeSize;
    long reads = 0;

    final BigLong2ShortHashMap hm; // for debug output

    public ReadsDispatcher(Source<Dna> reader, int workRangeSize, BigLong2ShortHashMap hmForMonitoring) {
        this.iterator = reader.iterator();
        this.workRangeSize = workRangeSize;
//        Tool.debug(logger, "Using " + workRangeSize + " reads as workRangeSize");
        this.hm = hmForMonitoring;
    }


    public synchronized List<Dna> getWorkRange() {
        List<Dna> list = new ArrayList<Dna>(workRangeSize);
        while ((list.size() < workRangeSize) && iterator.hasNext()) {
            list.add(iterator.next());
            ++reads;

            if (reads % 2500000 == 0) {
                Tool.debug(logger, "Processed " + NumUtils.groupDigits(reads) + " reads:");
                if (hm != null) {
                    Tool.debug(logger, "Total hm size = " + NumUtils.groupDigits(hm.size()) + ", " +
                            "size in hm.maps = {" + NumUtils.groupDigits(hm.maps[0].size()) + ", "
                            + NumUtils.groupDigits(hm.maps[1].size()) + ", "
                            + NumUtils.groupDigits(hm.maps[2].size()) + ", "
                            + NumUtils.groupDigits(hm.maps[3].size()) + ", ...}");
                }
                Tool.debug(logger, "Available memory (without running GC) = " + Misc.availableMemoryWithoutRunningGCAsString());
            }
        }
        return list.isEmpty() ? null : list;
    }
}
