import java.io.IOException;
import java.nio.file.FileSystems;
// import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.Paths;

public class TestPepIndex {

  public static void main(String[] args) throws IOException {
    Path path = Paths
        .get("input/test/20161021_ProAI_8nm_EDC_01.validated.intra.txt");
    PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:*.txt");
    System.out.println(path.getFileName());
    System.out.println(path.getNameCount());
    System.out.println(path.getName(1));
    System.out.println(path.getParent());
    System.out.println(matcher.matches(path.getFileName()));

    // Files.createDirectories(Paths.get("test1/test2/test3"));
  }
}
