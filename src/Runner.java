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


/*

        boolean printTools = (args.length > 0) && (args[0].equals("-ts") || args[0].equals("--tools"));
        if (printTools) {
            out.println("Available tools:");
            out.println();
            for (Tool t : tools) {
                out.println(fit(t.name, 30) + " " + fit(t.getClass().getName(), 40) + " " + t.description);
            }
            out.println();
            return;
        }

        boolean toolIsSet = (args.length > 0) && (args[0].equals("-t") || args[0].equals("--tool"));
        if (toolIsSet) {
            if (args.length < 2) {
                throw new RuntimeException("Tool name isn't specified");
            }

            Tool toolInst = null;
            for (Tool tool : tools) {
                if (tool.name.equals(args[1])) {
                    toolInst = tool;
                }
            }
            if (toolInst == null) {
                throw new RuntimeException("Tool name is incorrect");
            }

            args = Arrays.copyOfRange(args, 2, args.length);
            toolInst.mainImpl(args);
            return;
        }

        out.println("Runner: unknown parameters");
    }


    */

    public static void main(String[] args) {
        new Runner().run(args);
    }
}
