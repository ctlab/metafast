package algo;

import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import ru.ifmo.genetics.executors.BlockingThreadPoolExecutor;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import structures.Sequence;

import java.util.*;
import java.util.concurrent.ConcurrentLinkedDeque;

public class SequencesFinders {

    public static Deque<Sequence> thresholdStrategy(ArrayLong2IntHashMap hm,
                                                   int availableProcessors,
                                                   int freqThreshold,
                                                   int lenThreshold,
                                                   int k) throws InterruptedException {
        Deque<Sequence> ans = new ConcurrentLinkedDeque<Sequence>();
        LongOpenHashSet used = new LongOpenHashSet();

        BlockingThreadPoolExecutor executor = new BlockingThreadPoolExecutor(availableProcessors);

        for (int i = 0; i < hm.hm.length; ++i) {
            executor.blockingExecute(new
                    AddSequencesShiftingRightTask(hm, hm.hm[i], k, freqThreshold, lenThreshold, ans, used));
        }

//        System.out.println(executor.getTaskCount());
        executor.shutdownAndAwaitTermination();
        return ans;
    }
}
