import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * Find indexes for peptide 1 and peptide 2.
 * 
 * @author Zexi Jesse Zhuang
 */
public class PepIndex {
  /**
   * map protein Id to protein peptide sequence.
   */
  private static Map<String, String> proteins;

  private static void readProteins(File proteinFile) {
    proteins = new HashMap<>();
    try (Scanner scan = new Scanner(proteinFile)) {
      String firstLine;
      while (scan.hasNextLine() && (firstLine = scan.nextLine()) != null
          && firstLine.startsWith(">")) {
        int proteinIdStartIndex = firstLine.indexOf("|");
        if (proteinIdStartIndex == -1)
          System.out.println("Cannot find | in line starting with >");
        String proteinID = firstLine.substring(proteinIdStartIndex + 1,
            firstLine.indexOf("|", proteinIdStartIndex + 1));
        StringBuilder proteinPeptides = new StringBuilder();
        String line;
        while ((line = scan.nextLine()) != null && line.length() > 0)
          proteinPeptides.append(line);
        proteins.put(proteinID, proteinPeptides.toString());
      }

    } catch (FileNotFoundException e) {
      System.out.println("Cannot open protein file in input folder.");
    }
    System.out.println(proteins);
  }

  private static void findPeptidePositions(String peptideFile) {
    final String outputDirectory = "output/";
    final String delimiter = "\t";
    final String columns = delimiter + "peptide1" + delimiter + "position1"
        + delimiter + "peptide2" + delimiter + "position2";

    try (BufferedReader br = Files
        .newBufferedReader(Paths.get("input/" + peptideFile))) {
      try (BufferedWriter bw = Files.newBufferedWriter(
          Paths.get(outputDirectory + peptideFile), StandardCharsets.UTF_8,
          StandardOpenOption.CREATE, StandardOpenOption.WRITE)) {

        bw.write(br.readLine() + columns);
        bw.newLine();

        String line;
        while ((line = br.readLine()) != null) {
          int index = line.indexOf("-.");
          String peptide1 = line.substring(index + 2,
              (index = line.indexOf("(", index + 2)));

          int peptide1Pos = Integer
              .valueOf(line.substring(index + 1, line.indexOf(")", index)));
          index = line.indexOf("--");
          String peptide2 = line.substring(index + 2,
              (index = line.indexOf("(", index + 2)));
          int peptide2Pos = Integer
              .valueOf(line.substring(index + 1, line.indexOf(")", index)));

          index = line.indexOf("sp|");
          final String proteinId = line.substring(index + 3,
              line.indexOf("|", index + 3));
          final String proteinPeptides = proteins.get(proteinId);

          bw.write(line + delimiter + peptide1 + delimiter
              + (proteinPeptides.indexOf(peptide1) + peptide1Pos) + delimiter
              + peptide2 + delimiter
              + (proteinPeptides.indexOf(peptide2) + peptide2Pos));
          bw.newLine();
        }

      } catch (IOException e) {
        System.out.println("Cannot write output file.");
      }
    } catch (IOException e1) {
      System.out.println("Cannot read peptide file.");
    }
  }

  public static void main(String[] args) {
    // read configuration file
    final String inputDirectory = "input/";
    File configuration = new File(inputDirectory + "configuration.txt");
    try (Scanner scan = new Scanner(configuration)) {
      String proteinFile, peptideFile;
      while (scan.hasNextLine()) {
        proteinFile = scan.next();
        // read all proteins
        readProteins(new File(inputDirectory + proteinFile));

        // process experiment files, find indexes, and write output
        while (scan.hasNext()) {
          peptideFile = scan.next();
          findPeptidePositions(peptideFile);
        }
      }
    } catch (FileNotFoundException e) {
      System.out.println("Cannot open configuration file in input folder.");
    }

  }
}
