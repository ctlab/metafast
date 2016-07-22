package io;

import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.io.sources.Source;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.iterators.ProgressableIterator;
import ru.ifmo.genetics.utils.tool.Tool;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

public class BytesDispatcher {
    final Logger logger = Logger.getLogger("bytes-dispatcher");

    final InputStream is;
    public final int workRangeSize;
    long bytesRead = 0;

    final BigLong2ShortHashMap hm; // for debug output

    public BytesDispatcher(InputStream is, int workRangeSize, BigLong2ShortHashMap hmForMonitoring) {
        this.is = is;
        this.workRangeSize = workRangeSize;
//        Tool.debug(logger, "Using " + workRangeSize + " bytes as workRangeSize");
        hm = hmForMonitoring;
    }



    public byte[] getNewEmptyWorkRange() {
        return new byte[workRangeSize];
    }

    public synchronized int readWorkRange(byte[] range) {
        try {
            int read = 0;
            while (read < range.length) {
                int r = is.read(range, read, range.length - read);
                if (r == -1) {
                    break;
                }
                read += r;
            }

            bytesRead += read;
            /*
            if ((bytesRead & ((1 << 29) - 1)) == 0) { // 512 Mb
                Tool.debug(logger, "Processed " + (bytesRead >> 20) + " Mb of data:");
                if (hm != null) {
                    Tool.debug(logger, "Total hm size = " + NumUtils.groupDigits(hm.size()) + ", " +
                            "size in hm.maps = {" + NumUtils.groupDigits(hm.maps[0].size()) + ", "
                            + NumUtils.groupDigits(hm.maps[1].size()) + ", "
                            + NumUtils.groupDigits(hm.maps[2].size()) + ", "
                            + NumUtils.groupDigits(hm.maps[3].size()) + ", ...}");
                }
                Tool.debug(logger, "Available memory (without running GC) = " + Misc.availableMemoryWithoutRunningGCAsString());
            }
            */
            return read;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
