import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SAIS {

    public static void main(String[] args) throws IOException {

        // Default to alphabet size of 256, add 1 because of artificial shift of byte value by one for adding sentinel
        int alphabetSize = 257;

        byte[] inputBytes = Files.readAllBytes(Paths.get(args[0]));

        // Preprocessing; convert input to array of unsigned integers because Java bytes are signed, add one to all byte
        // values to ensure sentinel with byte value zero is smaller than all other characters appearing
        int[] input = new int[inputBytes.length + 1];
        for (int i = 0; i < inputBytes.length; i++) {
            input[i] = (((int) inputBytes[i]) & 0xff) + 1;
        }

        // Add sentinel with byte value zero at the very end
        input[inputBytes.length] = 0;


        ///////////////////
        // SUFFIX ARRAY
        ///////////////////

        // Get used memory and start time BEFORE SA construction
        long beforeUsedMem = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        long startTimeSA = System.nanoTime();

        // Construct suffix array with cyclic shift algorithm in O(n log n) time
        int[] sa = sais(input, alphabetSize);

        // Get used memory and start time AFTER SA construction
        long endTimeSA = System.nanoTime();
        long afterUsedMem = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();

        // Calculate used time and memory of SA construction
        long timeUsedSA = (endTimeSA - startTimeSA) / 1000000;
        long memUsed = (afterUsedMem - beforeUsedMem) / 1000000;


        ///////////////////
        // LCP NAIVE
        ///////////////////
        long startTimeLCPNaive = System.nanoTime();
        int[] naiveLCP = naiveLCPArray(input, sa);
        long endTimeLCPNaive = System.nanoTime();

        long timeUsedLCPNaive = (endTimeLCPNaive - startTimeLCPNaive) / 1000000;


        ///////////////////
        // LCP KASAI
        ///////////////////
        long startTimeLCPKasai = System.nanoTime();
        int[] kasaiLCP = kasaiLCPArray(input, sa);
        long endTimeLCPKasai = System.nanoTime();

        long timeUsedLCPKasai = (endTimeLCPKasai - startTimeLCPKasai) / 1000000;


        ///////////////////
        // LCP PHI
        ///////////////////
        long startTimeLCPPhi = System.nanoTime();
        int[] phiLCP = phiLCPArray(input, sa);
        long endTimeLCPPhi = System.nanoTime();

        long timeUsedLCPPhi = (endTimeLCPPhi - startTimeLCPPhi) / 1000000;


        // Final program output print
        System.out.printf("RESULT name=TobiasSchiebel sa_construction_time=%d sa_construction_memory=%d lcp_naive_construction_time=%d lcp_kasai_construction_time=%d lcp_phi_construction_time=%d", timeUsedSA, memUsed, timeUsedLCPNaive, timeUsedLCPKasai, timeUsedLCPPhi);
    }


    /**
     * Computes the suffix array of a given text using the SA-IS algorithm.
     *
     * @param text         The given text.
     * @param alphabetSize The size of the underlying alphabet of the text.
     * @return the resulting suffix array.
     */
    private static int[] sais(int[] text, int alphabetSize) {

        // Classify every suffix as L (0), S (1) or S*-type (2)
        int[] classes = classification(text);

        // Calculate the sizes of the respective buckets for each character in the alphabet to store buckets more
        // efficiently as both the head and tail of a bucket can be computed from the given sizes
        int[] bucketSizes = calcBucketSizes(text, alphabetSize);

        // Naively put LMS suffixes (S*-type) in suffix array at end of buckets in text order
        int[] suffixArray = insertLMSNaive(text, bucketSizes, classes);

        // Assign all other suffixes according to naively placed LMS suffixes
        induceL(text, suffixArray, bucketSizes, classes);
        induceS(text, suffixArray, bucketSizes, classes);

        // Construct a reduced text of meta symbols from positions of LMS substrings in the suffix array because LMS
        // suffixes are not necessarily sorted correctly yet
        int[][] reduceResult = reduceSuffixArray(text, suffixArray, classes);

        // Expand reduction results
        int[] reducedText = reduceResult[0];
        int reducedAlphabetSize = reduceResult[1][0];
        int[] reducedOffsets = reduceResult[2];

        // Construct the suffix array for the reduced text, potentially go into recursion if reduced text is too complex
        int[] reducedSuffixArray = getReducedSuffixArray(reducedText, reducedAlphabetSize);

        // Put LMS substrings to their correct positions with remapping from suffix array of reduced text
        suffixArray = insertLMSExact(text, bucketSizes, reducedSuffixArray, reducedOffsets);

        // Assign all other suffixes according to accurately placed LMS suffixes
        induceL(text, suffixArray, bucketSizes, classes);
        induceS(text, suffixArray, bucketSizes, classes);

        return suffixArray;
    }


    /**
     * Assigns classes S (smaller) or L (larger) to each character in input text dependent on whether their successor
     * is lexicographically larger.
     * The class S becomes S* (leftmost S suffix) whenever the predecessor is lexicographically larger.
     *
     * @param text The input text as unsigned byte codes with length n.
     * @return int[] The array of class representations as L=0, S=1, S*=2 of length n.
     */
    private static int[] classification(int[] text) {
        int n = text.length;

        // If there are no more characters left to classify classification is finished immediately
        if (n == 0) {
            return new int[0];
        }

        int[] classes = new int[n];

        // Assign S and L types from right to left
        for (int i = n - 1; i >= 0; i--) {
            // S type
            if (i == n - 1 || text[i] < text[i + 1]) {
                classes[i] = 1;

            }
            // L type
            else if (text[i] > text[i + 1]) {
                classes[i] = 0;
                // S* type when previous character had L type
                if (classes[i + 1] == 1) {
                    classes[i + 1] = 2;
                }

            }
            // Subsequent characters
            else {
                classes[i] = classes[i + 1];
            }
        }
        return classes;
    }


    /**
     * Determines bucket sizes for the given text.
     *
     * @param text         The given text.
     * @param alphabetSize The size of the underlying alphabet of the text.
     * @return how large the bucket for the char at the given index is aka how often the char appears in the text
     */
    private static int[] calcBucketSizes(int[] text, int alphabetSize) {
        int[] bucketSizes = new int[alphabetSize];
        Arrays.fill(bucketSizes, 0);

        // Count how often the char appears in the text
        for (int c : text) {
            bucketSizes[c]++;
        }

        return bucketSizes;
    }


    /**
     * Determines where the bucket for each character starts given an array of bucket sizes.
     *
     * @param bucketSizes The given bucket sizes for a text.
     * @return where the bucket of a char starts.
     */
    private static int[] getBucketHeads(int[] bucketSizes) {
        int offset = 0;
        int[] bucketHeads = new int[bucketSizes.length];

        for (int i = 0; i < bucketSizes.length; i++) {
            bucketHeads[i] = offset;
            offset += bucketSizes[i];
        }

        return bucketHeads;
    }


    /**
     * Determines where the bucket for each character ends given an array of bucket sizes.
     *
     * @param bucketSizes The given bucket sizes for a text.
     * @return where the bucket of a char ends.
     */
    private static int[] getBucketTails(int[] bucketSizes) {
        int offset = 0;
        int[] bucketTails = new int[bucketSizes.length];

        for (int i = 0; i < bucketSizes.length; i++) {
            offset += bucketSizes[i];
            bucketTails[i] = offset - 1;
        }

        return bucketTails;
    }


    /**
     * Inserts LMS (S*) suffixes naively into the suffix array at the end of their character's bucket in text order.
     *
     * @param text        The given text.
     * @param bucketSizes The bucket sizes of the given text.
     * @param classes     The classes of the given text.
     * @return the suffix array initialized with -1 and the S*-suffixes inserted naively.
     */
    private static int[] insertLMSNaive(int[] text, int[] bucketSizes, int[] classes) {

        int n = text.length;

        // Initialize suffix array with -1 as placeholder
        int[] naiveSA = new int[n];
        Arrays.fill(naiveSA, -1);

        int[] bucketTails = getBucketTails(bucketSizes);

        for (int i = 0; i < n; i++) {
            // Only consider S*-type suffixes
            if (classes[i] != 2) {
                continue;
            }

            // Assign to end of current char bucket
            int currentChar = text[i];
            naiveSA[bucketTails[currentChar]] = i;

            // Decrease bucket size at end of current char bucket
            bucketTails[currentChar]--;
        }

        return naiveSA;
    }


    /**
     * Scan suffix array from left to right and assign the L-type suffixes.
     *
     * @param text        The given text.
     * @param suffixArray The given suffix array (initialized with -1 and LMS suffixes in place).
     * @param bucketSizes The bucket sizes of the given text.
     * @param classes     The classes of the given text.
     */
    private static void induceL(int[] text, int[] suffixArray, int[] bucketSizes, int[] classes) {
        int[] bucketHeads = getBucketHeads(bucketSizes);

        for (int i = 0; i < suffixArray.length; i++) {
            // Skip placeholder values
            if (suffixArray[i] == -1) {
                continue;
            }

            // Suffix left of the entry linked in suffix array
            int left = suffixArray[i] - 1;

            // Make sure index is in bounds and left suffix is L-type
            if (left >= 0 && classes[left] == 0) {
                // Put left suffix to start of its bucket and decrease bucket size at start of left char bucket
                int leftChar = text[left];
                suffixArray[bucketHeads[leftChar]] = left;
                bucketHeads[leftChar]++;
            }
        }
    }


    /**
     * Scan suffix array from right to left and assign the S-type suffixes.
     *
     * @param text        The given text.
     * @param suffixArray The given suffix array (initialized with -1 and LMS suffixes in place).
     * @param bucketSizes The bucket sizes of the given text.
     * @param classes     The classes of the given text.
     */
    private static void induceS(int[] text, int[] suffixArray, int[] bucketSizes, int[] classes) {
        int[] bucketTails = getBucketTails(bucketSizes);

        for (int i = suffixArray.length - 1; i >= 0; i--) {
            // Suffix left of the entry linked in suffix array
            int left = suffixArray[i] - 1;

            // Make sure index is in bounds and left suffix is S-type
            if (left >= 0 && classes[left] != 0) {
                // Put left suffix to start of its bucket and decrease bucket size at start of left char bucket
                int leftChar = text[left];
                suffixArray[bucketTails[leftChar]] = left;
                bucketTails[leftChar]--;
            }
        }
    }


    /**
     * Checks whether two LMS substrings at offsets offsetA and offsetB are equal
     *
     * @param text    The text to check.
     * @param classes The classification of the text.
     * @param offsetA The offset where the first LMS substring starts.
     * @param offsetB The offset where the second LMS substring starts.
     * @return true if the LMS substrings at the given offsets are equal.
     */
    private static boolean lmsEquals(int[] text, int[] classes, int offsetA, int offsetB) {
        if (offsetA == text.length || offsetB == text.length) {
            return false;
        }

        int i = 0;
        while (true) {
            // Check whether the current characters looked at are S*-type (sentinel at the end is also always S*,
            // therefore we will always end up in´one of the first two cases and the loop will terminate when one
            // substring reaches the end of the text)
            boolean aIsLMS = classes[offsetA + i] == 2;
            boolean bIsLMS = classes[offsetB + i] == 2;

            // Reached end of LMS substrings
            if (i > 0 && aIsLMS && bIsLMS) {
                return true;
            }

            // Different classes
            if (aIsLMS != bIsLMS) {
                return false;
            }

            // Different characters in text
            if (text[i + offsetA] != text[i + offsetB]) {
                return false;
            }

            i++;
        }
    }


    /**
     * Constructs a reduced text of meta symbols from positions of LMS substrings in the suffix array.
     *
     * @param text        The given text.
     * @param suffixArray The suffix array with naively placed LMS substrings and after inducing L & S.
     * @param classes     The classes of the given text.
     * @return the result of the reduction consisting of [reducedText, reducedAlphabetSize, reducedOffsets].
     */
    private static int[][] reduceSuffixArray(int[] text, int[] suffixArray, int[] classes) {
        int n = text.length;

        // Array to store the new meta symbols (names) of the LMS substrings
        int[] names = new int[n];
        Arrays.fill(names, -1);

        // The incremental names we assign to LMS substrings starting from zero
        int nameIndex = 0;

        // Offset of the last LMS suffix that was checked
        int lastLMSOffset = -1;

        for (int i = 0; i < suffixArray.length; i++) {
            int suffixOffset = suffixArray[i];

            // Only LMS suffixes are relevant
            if (classes[suffixOffset] != 2) {
                continue;
            }

            // Assign new name when LMS substrings are different
            if (lastLMSOffset != -1 && !lmsEquals(text, classes, lastLMSOffset, suffixOffset)) {
                nameIndex++;
            }

            // Update last LMS suffix to current
            lastLMSOffset = suffixOffset;

            // Save name to names array
            names[suffixOffset] = nameIndex;
        }

        List<Integer> reducedOffsets = new ArrayList<>();
        List<Integer> reducedText = new ArrayList<>();

        for (int i = 0; i < names.length; i++) {
            // Ignore placeholder names
            if (names[i] == -1) {
                continue;
            }

            // Store which original LMS substring is represented in the reduced text for remapping later
            reducedOffsets.add(i);
            // Put together meta symbols in text order
            reducedText.add(names[i]);
        }
        return new int[][]{reducedText.stream().mapToInt(Integer::intValue).toArray(), {nameIndex + 1}, reducedOffsets.stream().mapToInt(Integer::intValue).toArray()};
    }


    /**
     * Constructs the suffix array of a given reduced text.
     *
     * @param reducedText         The reduced text consisting of meta symbols.
     * @param reducedAlphabetSize The alphabet size of the reduced text.
     * @return the suffix array of the reduced text.
     */
    private static int[] getReducedSuffixArray(int[] reducedText, int reducedAlphabetSize) {
        int[] reducedSuffixArray;

        // Every character in the reduced text appears only once
        if (reducedAlphabetSize == reducedText.length) {
            reducedSuffixArray = new int[reducedText.length + 1];
            Arrays.fill(reducedSuffixArray, -1);

            // Add sentinel at beginning
            reducedSuffixArray[0] = reducedText.length;

            // Bucket sort, the incremental naming of the meta symbols ensures the correct buckets are assigned.
            // Bucket heads / tails not required since every symbol appears only once.
            for (int i = 0; i < reducedText.length; i++) {
                int metaSymbol = reducedText[i];
                reducedSuffixArray[metaSymbol + 1] = i;
            }
        }
        // Reduced text is still too complex, thus recursion is required
        else {
            reducedSuffixArray = sais(reducedText, reducedAlphabetSize);
        }

        return reducedSuffixArray;
    }


    /**
     * Inserts LMS (S*) suffixes accurately into the suffix array at the end of their character's buckets in the correct
     * order by remapping from the suffix array of the reduced text.
     *
     * @param text               The original given text.
     * @param bucketSizes        The bucket sizes of the given text.
     * @param reducedSuffixArray The suffix array of the reduced text.
     * @param reducedOffsets     The offsets of the reduced text required for remapping to the original text.
     * @return the suffix array of the text with the LMS suffixes at their accurate positions.
     */
    private static int[] insertLMSExact(int[] text, int[] bucketSizes, int[] reducedSuffixArray, int[] reducedOffsets) {
        int[] suffixArray = new int[text.length];
        Arrays.fill(suffixArray, -1);

        // The LMS suffixes are added at the end of the buckets, thus iterate over reduced SA in reverse
        int[] bucketTails = getBucketTails(bucketSizes);
        for (int i = reducedSuffixArray.length - 1; i >= 1; i--) {
            // Remapping to char in text from reduced SA with offsets
            int textIndex = reducedOffsets[reducedSuffixArray[i]];
            int charInText = text[textIndex];

            // Assign to end of char bucket and decrease bucket size at end of bucket
            suffixArray[bucketTails[charInText]] = textIndex;
            bucketTails[charInText]--;
        }

        suffixArray[0] = text.length - 1;

        return suffixArray;
    }


    /**
     * Computes the LCP-Array naively in quadratic time.
     *
     * @param text The input text.
     * @param sa   The suffix array of the text.
     * @return the LCP-Array of the text.
     */
    private static int[] naiveLCPArray(int[] text, int[] sa) {
        // Initialize LCP-Array
        int[] naiveLCP = new int[text.length];
        naiveLCP[0] = 0;

        // Compare all suffixes to previous suffix
        for (int i = 1; i < text.length; i++) {
            naiveLCP[i] = compareLCP(text, sa[i - 1], sa[i]);
        }
        return naiveLCP;
    }


    /**
     * Calculates the length of the longest common prefix (LCP) of two suffixes starting from offsetA and offsetB for
     * naive LCP-Array computation.
     *
     * @param text    The input text.
     * @param offsetA The offset of the first suffix.
     * @param offsetB The offset of the second suffix.
     * @return the length of the longest common prefix of the suffixes at offsetA and offsetB.
     */
    private static int compareLCP(int[] text, int offsetA, int offsetB) {
        int i = 0;
        while (text[i + offsetA] == text[i + offsetB]) {
            i++;
        }
        return i;
    }


    /**
     * Computes the LCP-Array in linear time with the algorithm from Kasai et al. using the inverse suffix array.
     *
     * @param text The input text.
     * @param sa   The suffix array of the text.
     * @return the LCP-Array of the text.
     */
    private static int[] kasaiLCPArray(int[] text, int[] sa) {
        // Initialize LCP-Array
        int[] kasaiLCP = new int[text.length];
        kasaiLCP[0] = 0;

        // Compute inverse suffix array
        int[] inverseSA = getInverseSuffixArray(sa);

        // Iterate over inverse suffix array
        int lcpLength = 0;
        for (int i = 0; i < text.length; i++) {

            // Ensure index is in bounds
            if (inverseSA[i] == 0) {
                continue;
            }

            // Offset of suffix to the left of current suffix
            int j = sa[inverseSA[i] - 1];

            // Compare suffixes at offsets i and j and count common characters
            while (text[i + lcpLength] == text[j + lcpLength]) {
                lcpLength++;
            }

            // Assign LCP-Array value and reset LCP length counter
            kasaiLCP[inverseSA[i]] = lcpLength;
            lcpLength = Math.max(0, lcpLength - 1);
        }
        return kasaiLCP;
    }


    /**
     * Computes the inverse suffix array containing the inverse permutations (where a suffix is in the suffix array).
     *
     * @param sa The suffix array.
     * @return the inverse suffix array.
     */
    private static int[] getInverseSuffixArray(int[] sa) {
        int[] inverseSA = new int[sa.length];

        for (int i = 0; i < sa.length; i++) {
            inverseSA[sa[i]] = i;
        }

        return inverseSA;
    }


    /**
     * Computes the LCP-Array in linear time with the algorithm from Kärkkäinen et al. using the Phi-Array.
     *
     * @param text The input text.
     * @param sa   The suffix array of the text.
     * @return the LCP-Array of the text.
     */
    private static int[] phiLCPArray(int[] text, int[] sa) {
        // Initialize LCP-Array
        int[] phiLCP = new int[text.length];

        // Compute Phi-Array
        int[] phi = getPhiArray(sa);

        // Iterate over Phi-Array (text order)
        int lcpLength = 0;
        for (int i = 0; i < text.length; i++) {

            // Get offset of suffix for comparison from Phi-Array
            int j = phi[i];

            // Compare suffixes at offsets i and j and count common characters
            while (text[i + lcpLength] == text[j + lcpLength]) {
                lcpLength++;
            }

            // Assign LCP value at current index and reset LCP length counter
            phi[i] = lcpLength;
            lcpLength = Math.max(0, lcpLength - 1);
        }

        // Reorder to obtain suffix array order
        for (int i = 0; i < text.length; i++) {
            phiLCP[i] = phi[sa[i]];
        }

        return phiLCP;
    }


    /**
     * Computes the Phi-Array containing the suffix required for comparison when building the LCP-Array.
     *
     * @param sa The suffix array.
     * @return the inverse suffix array.
     */
    private static int[] getPhiArray(int[] sa) {
        int[] phi = new int[sa.length];

        for (int i = 1; i < sa.length; i++) {
            phi[sa[i]] = sa[i - 1];
        }

        return phi;
    }
}
