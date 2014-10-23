public class Runner extends ru.ifmo.genetics.Runner {
    @Override
    protected void printHeader() {
        out.println("Fast metagenome analysis toolkit, version " + getVersion());
        out.println();
    }

    @Override
    protected void printHelp() {
        out.println("Usage: java [<JVM options>] -jar <path-to-jar>/metafast.jar [<Launch options>] [<Tool parameters>]");
        out.println("");
        out.println("This toolkit allows you to run different tools from it.");
        out.println("To see available tools:              java -jar <path-to-jar>/metafast.jar -ts");
        out.println("To see help for selected tool:       java -jar <path-to-jar>/metafast.jar -t <tool-name>");
        out.println("To run selected tool:                java -jar <path-to-jar>/metafast.jar -t <tool-name> <tool-parameters>");
        out.println("");
    }

    @Override
    protected void checkOptionsBeforeRunning(String[] args) {
        boolean printVersion = containsOption(args, "--version");
        boolean printHelp = (args.length == 0) ||
                (containsOption(args, "-h", "--help") && !containsOption(args, "-t", "--tool"));

        if (printVersion || printHelp) {
            printHeader();
        }
        if (printHelp) {
            printHelp();
        }
        if (printVersion || printHelp) {
            System.exit(0);
        }
    }



    public static void main(String[] args) {
        new Runner().run(args);
    }
}
