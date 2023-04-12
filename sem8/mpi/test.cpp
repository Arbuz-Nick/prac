#include <iostream>
int main(int argc, char const *argv[])
{
    int rank, size;
    sscanf(argv[1], "%d", &rank);
    sscanf(argv[2], "%d", &size);
    int k = -1;
    for (int i = 0; i < 9; i++)
    {
        if (i % (size - 1) == 0)
        {
            k++;
        } 
        std::cout << i << ' ' << (rank + 1 + (i % (size - 1))) % size << std::endl;
    }

    return 0;
}
