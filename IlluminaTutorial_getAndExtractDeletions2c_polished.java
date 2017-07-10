//package net.sf.samtools.example;

import net.sf.samtools.*;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.Scanner;
import java.io.File;

//filter out unmapped reads before using

public class IlluminaTutorial_getAndExtractDeletions2c
{
	static private final int MAX_ACCEPTABLE_INSERT = 875; //419 + 152*3 --> the median insert size plus three standard deviations, found using InferredInsertSizeMetrics
	static private String inputFile = "NA12878.hs37d5.bam";
	static private String outputFile = "NA12878.hs37d5_deletionsV2.bam";
	static private final int INTERMEDIATE_READS_THRESHOLD = 0; //how many reads are allowable within an identified deletion

	public static void main(String args[])
	{
		File inputSamOrBamFile = new File(inputFile);
		File outputSamOrBamFile = new File(outputFile);


		final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);

		int currentInsertLength = 0;
		int nextInsertLength = 0;
		String currentChromosome = "1";
		String nextChromosome = "1";
		int currentStart = 0;
		int nextAlignmentStart = 0;
		int nextCigarLength = 0;
		int nextStart = 0;
		int currentEnd = 0;
		int nextEnd = 0;
		int currentDeletionAlignments = 0;
		boolean potentialDeletion = false; //if you are currently investigating a potential deletion


		System.out.println("Chromosome  |  Start  |  End  |  Length  |  Instances");

		for (final SAMRecord samRecord : inputSam)
		{
			nextChromosome = samRecord.getReferenceName();
			nextInsertLength = samRecord.getInferredInsertSize();
			nextAlignmentStart = samRecord.getAlignmentStart();
			nextCigarLength = samRecord.getCigarLength();
			nextStart = nextAlignmentStart + nextCigarLength;
			nextEnd = samRecord.getMateAlignmentStart();
			
			if(nextInsertLength > MAX_ACCEPTABLE_INSERT && nextEnd > nextAlignmentStart)
			{//the insert size of the alignment and its mate is greater than the acceptable length, indicating a potential deletion
				potentialDeletion = true;
				if(nextChromosome == currentChromosome && currentInsertLength/nextInsertLength > 0.8 && currentInsertLength/nextInsertLength < 1.25 && nextStart < currentEnd)
				{//if the new potential deletion is the same as the one that is currently being monitored (same chromosome, similar insert lengths, and same range on reference) - assuming alignments are coordinate-sorted; otherwise, must check if nextStart > currentStart and currentEnd < nextEnd **think about when paired reads are flipped
					currentStart = nextStart;
					currentDeletionAlignments++;
					if(currentEnd > nextEnd)
						currentEnd = nextEnd;
				}
				else
				{//if the new potential deletion is a new one
					if(currentDeletionAlignments > 0)
					{//if this also marks the end of the previous deletion
						if((currentStart + 1000) < currentEnd && currentEnd <= nextAlignmentStart)
							System.out.println(currentChromosome + "\t" + currentStart + "\t" + currentEnd + "\t" + (currentEnd - currentStart + 1) + "\t" + currentDeletionAlignments);
						else if((currentStart + 1000) < nextAlignmentStart && nextStart < currentEnd)
							System.out.println(currentChromosome + "\t" + currentStart + "\t" + nextAlignmentStart + "\t" + (nextAlignmentStart - currentStart + 1) + "\t" + currentDeletionAlignments);
					}
					currentChromosome = nextChromosome;
					currentStart = nextStart;
					currentEnd = nextEnd;
					currentInsertLength = nextInsertLength;
					currentDeletionAlignments = 1;
				}
			}

			else if (nextEnd > nextAlignmentStart && nextInsertLength >= 0) //looking at reads aligned in the correct library orientation only
			{
				if(potentialDeletion == true && nextAlignmentStart > (currentStart + 1000))
				{//jump in alignment starts marks the confirmed location of deletion while a potential deletion is being investigated
					if (samRecord.getAlignmentStart() > currentEnd)
						System.out.println(currentChromosome + "\t" + currentStart + "\t" + currentEnd + "\t" + (currentEnd - currentStart + 1) + "\t" + currentDeletionAlignments);
					else
						System.out.println(currentChromosome + "\t" + currentStart + "\t" + nextAlignmentStart + "\t" + (nextAlignmentStart - currentStart + 1) + "\t" + currentDeletionAlignments);
					
					currentStart = 0;
					currentEnd = 0;
					currentInsertLength = 0;
					currentDeletionAlignments = 0;
					potentialDeletion = false;
				}
				else if(potentialDeletion == true && nextAlignmentStart > currentEnd)
				{//the current deletion is confirmed by the start of the next alignment exceeding the current predicted end of the deletion
					System.out.println(currentChromosome + "\t" + currentStart + "\t" + currentEnd + "\t" + (currentEnd - currentStart + 1) + "\t" + currentDeletionAlignments);
					currentStart = 0;
					currentEnd = 0;
					currentInsertLength = 0;
					currentDeletionAlignments = 0;
					potentialDeletion = false;
				}


				else if (potentialDeletion == true)
				{//the start and end predictions of the deletion are adjusted based on what is read
					if (currentStart < nextStart)
						currentStart = nextStart;
					if (currentEnd > nextEnd)
						currentEnd = nextEnd;
				}
			}

			if ((currentStart + 1000) > currentEnd && potentialDeletion == true)
			{//the flagged potential deletion is too short to be considered an SV, so values are reset
				potentialDeletion = false;
				currentStart = 0;
				currentEnd = 0;
				currentInsertLength = 0;
				currentDeletionAlignments = 0;
			}

		}

		inputSam.close();

	}
}
