#include "transcript.h"
#include "utils.h"
#include "string.h"
#include <stdio.h>

#define BUFFER_SIZE 10000
char buffer[BUFFER_SIZE];

struct ref_info_t read_ref_info(char *filename)
{
        struct ref_info_t ref;

        ref.name = NULL;
        ref.len = NULL;

        FILE *f = fopen(filename, "r");

        int count = 0;
        while (1) {
                if (!fgets(buffer, BUFFER_SIZE, f))
                        break;
                
                if (buffer[0] == '>') { 
                        ++count;
                        ref.name = realloc(ref.name, count * sizeof(*ref.name));
                        ref.len = realloc(ref.len, count * sizeof(*ref.len));

                        ref.name[count - 1] = strdup(buffer + 1);
                        ref.len[count - 1] = 0;
                }

                int len = strlen(buffer);

                if (buffer[len-1] == '\n')
                        --len;

                ref.len[count-1] += len;
        }

        fclose(f);

        ref.nref = count;
        return ref;
}


void release_ref_info(struct ref_info_t *ref)
{
        int i;
        for (i = 0; i < ref->nref; ++i)
                free(ref->name[i]);
        
        free(ref->name);
        free(ref->len);
}