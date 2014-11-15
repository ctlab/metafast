package io;

import java.util.concurrent.CountDownLatch;

public abstract class BytesWorker implements Runnable {

    private BytesDispatcher dispatcher = null;
    private CountDownLatch latch = null;

    boolean interrupted = false;


    void setDispatcher(BytesDispatcher dispatcher) {
        this.dispatcher = dispatcher;
    }
    void setLatch(CountDownLatch latch) {
        this.latch = latch;
    }



    public abstract void process(byte[] range, int len);

    @Override
    public void run() {
        if (dispatcher == null || latch == null) {
            throw new RuntimeException("Not full initialization!");
        }
        while (!interrupted) {
            byte[] range = dispatcher.getNewEmptyWorkRange();
            int r = dispatcher.readWorkRange(range);
            if (r <= 0) {
                break;
            }
            process(range, r);
        }
        latch.countDown();
    }

    public void interrupt() {
        interrupted = true;
    }
}
