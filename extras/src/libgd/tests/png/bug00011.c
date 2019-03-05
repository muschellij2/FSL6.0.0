/* $Id: bug00011.c,v 1.1.1.1 2014/02/06 18:14:22 duncan Exp $ */
#include "gd.h"
#include <stdio.h>
#include <stdlib.h>
#include "gdtest.h"

int main()
{
	gdImagePtr im;
	FILE *fp;
	char path[2048];

	sprintf(path, "%s/png/emptyfile", GDTEST_TOP_DIR);
	fp = fopen(path, "rb");
	if (!fp) {
		fprintf(stderr, "failed, cannot open file: %s\n", path);
		return 1;
	}
	im = gdImageCreateFromPng(fp);
	fclose(fp);

	if (!im) {
		return 0;
	} else {
		return 1;
	}
}
