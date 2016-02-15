/*-----------------------------------------------------------------------
|						     Count the lines         					|
-----------------------------------------------------------------------*/
int linesoffilex(FILE *estream)
{
	
	fseek(estream, 0, SEEK_SET);
	printf("step 1: Counting the lines of the file\n");
	
	char tx[12], ty[12], tz[12], tin[12];
	
	char* forpop = "%11s %11s %11s %11s\n";
	
	int j=0;
	
	while (fscanf (estream, forpop,
				   tx, ty, tz, tin) == 4)
	{
		j++;
	}
	
	printf("%d\n", j);
	return j;
}