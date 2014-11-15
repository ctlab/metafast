package algo;

import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import ru.ifmo.genetics.executors.BlockingThreadPoolExecutor;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import structures.Sequence;

import java.util.*;
import java.util.concurrent.ConcurrentLinkedDeque;

public class SequencesFinders {

    public static Deque<Sequence> thresholdStrategy(BigLong2ShortHashMap hm,
                                                   int availableProcessors,
                                                   int freqThreshold,
                                                   int lenThreshold,
                                                   int k) throws InterruptedException {
        Deque<Sequence> ans = new ConcurrentLinkedDeque<Sequence>();
        LongOpenHashSet used = new LongOpenHashSet();

        BlockingThreadPoolExecutor executor = new BlockingThreadPoolExecutor(availableProcessors);

        for (int i = 0; i < hm.maps.length; ++i) {
            executor.blockingExecute(new
                    AddSequencesShiftingRightTask(hm, hm.maps[i], k, freqThreshold, lenThreshold, ans, used));
        }

//        System.out.println(executor.getTaskCount());
        executor.shutdownAndAwaitTermination();
        return ans;
    }
}
