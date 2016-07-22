import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.StringParameterBuilder;
import ru.ifmo.genetics.utils.tool.parameters.ParameterDescription;
import tools.DistanceMatrixBuilderMain;

public class Runner extends ru.ifmo.genetics.Runner {

    public static final ParameterDescription<String> toolParameterDescr = new StringParameterBuilder("tool")
            .withShortOpt("t")
            .withDescription("set certain tool to run")
            .withDefaultComment(DistanceMatrixBuilderMain.NAME)
            .create();
    public static final ParameterDescription<String> memoryParameterDescr = new StringParameterBuilder("memory")
            .important()
            .withShortOpt("m")
            .withDescription("memory to use (for example: 1500M, 4G, etc.)")
            .withDefaultComment("90% of free memory (currently " + Misc.availableMemoryAsString() + ")")
            .create();

    static {
        ru.ifmo.genetics.Runner.toolParameter.replaceDescription(toolParameterDescr);
        Tool.launchOptions.remove(ru.ifmo.genetics.Runner.guiParameter);
        ru.ifmo.genetics.Runner.memoryParameter.replaceDescription(memoryParameterDescr);
    }

    Runner() {
        super();
        defaultTool = findTool(DistanceMatrixBuilderMain.class);
    }


    @Override
    protected void printHeader() {
        out.println("Fast metagenome analysis toolkit, version " + getVersion());
        out.println();
    }

    @Override
    protected void printFirstHelp() {
        out.println("Usage:     metafast [<Launch options>] [<Input parameters>]");
        printLater.add("To see full documentation visit https://github.com/ctlab/metafast/wiki");
    }

    @Override
    protected boolean checkOptionsBeforeRunning(String[] args) {
        boolean printVersion = containsOption(args, "--version");
        boolean printFirstHelp = (args.length == 0) ||
                (containsOption(args, getOptKeys(Tool.helpParameter)) && !containsOption(args, getOptKeys(toolParameter)));

        if (printVersion) {
            printHeader();
            return true;
        }
        if (printFirstHelp) {
            printHeader();
            printFirstHelp();
            return false;
        }
        return false;
    }



    public static void main(String[] args) {
        new Runner().run(args);
    }
}
