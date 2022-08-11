#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main()
{
    char hostbuffer[256];
    int hostname;

    // To retrieve hostname
    hostname = gethostname(hostbuffer, sizeof(hostbuffer));

    printf("Hostname: %s\n", hostbuffer);

    return 0;
}
