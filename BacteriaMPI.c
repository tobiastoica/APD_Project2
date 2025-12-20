
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <string.h>

// Global variables for command line arguments
int MAXITER;
int NTHREADS;

//#define DEBUG

pthread_barrier_t barrier;
pthread_barrier_t barrier2;

pthread_t *threads;  // Dynamic allocation since NTHREADS is now a variable

char *grid;
char *new_grid;
char *groundtruth;
int rows, cols;
 
int read_from_file(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("Error opening file: %s\n", filename);
        return 0;
    }

    // Read dimensions into global variables
    if (fscanf(file, "%d %d", &rows, &cols) != 2)
    {
        printf("Error reading dimensions from file\n");
        fclose(file);
        return 0;
    }

    // Allocate memory for grids
    grid = (char *)malloc(rows * cols * sizeof(char));
    if (!grid)
    {
        printf("Memory allocation error for grid\n");
        fclose(file);
        return 0;
    }
    new_grid = (char *)malloc(rows * cols * sizeof(char));
    if (!new_grid)
    {
        printf("Memory allocation error for new grid\n");
        free(grid);
        fclose(file);
        return 0;
    }

    // Read grid configuration
    // Skip the newline after dimensions
    fgetc(file);
    
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            char ch = fgetc(file);
            // Skip newlines
            while (ch == '\n' || ch == '\r')
            {
                ch = fgetc(file);
            }
            
            if (ch == 'X' || ch == 'x')
            {
                grid[i * cols + j] = 1; // Bacteria present
            }
            else if (ch == '.')
            {
                grid[i * cols + j] = 0; // Empty
            }
            else if (ch == EOF)
            {
                printf("Unexpected end of file\n");
                free(grid);
                free(new_grid);
                fclose(file);
                return 0;
            }
        }
    }

    fclose(file);
    return 1;
}
 
void write_grid(const char *output_filename)
{
    FILE *file = fopen(output_filename, "w");
    if (!file)
    {
        printf("Error opening output file: %s\n", output_filename);
        return;
    }

    // Write dimensions
    fprintf(file, "%d %d\n", rows, cols);

    // Write grid
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (grid[i * cols + j] == 1)
            {
                fprintf(file, "X");
            }
            else
            {
                fprintf(file, ".");
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);
    printf("Grid saved to %s\n", output_filename);
}

void print_grid(void)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (grid[i * cols + j] == 1)
            {
                printf("X");
            }
            else
            {
                printf(".");
            }
        }
        printf("\n");
    }
}
 
 int equal_groundtruth(void)
 {
     for (int i = 0; i < rows; i++)
         for (int j = 0; j < cols; j++)
         {
             if (grid[i * cols + j] != groundtruth[i * cols + j])
                 return 0;
         }
     return 1;
 }
 
void save_groundtruth(void)
{
    groundtruth = (char *)malloc(rows * cols * sizeof(char));
    if (!groundtruth)
    {
        printf("Memory allocation error for groundtruth result\n");
        exit(1);
    }

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
        {
            groundtruth[i * cols + j] = grid[i * cols + j];
        }
}
 
void swap_ptr(char **p1, char **p2)
{
    char *tmp = *p1;
    *p1 = *p2;
    *p2 = tmp;
}
 
int number_of_neighbors(int i, int j)
{
    int count = 0;
    for (int di = -1; di <= 1; di++)
    {
        for (int dj = -1; dj <= 1; dj++)
        {
            if (di == 0 && dj == 0)
                continue;
            
            int ni = i + di;
            int nj = j + dj;
            
            if (ni >= 0 && ni < rows && nj >= 0 && nj < cols)
            {
                count += grid[ni * cols + nj];
            }
        }
    }
    return count;
}

void serial_bacteria()
{
    int i, j, time, neighbors;

    for (time = 0; time < MAXITER; time++)
    {
    #ifdef DEBUG
       printf("\nGeneration %d \n", time);
       print_grid();
    #endif
        // Process all cells including boundaries
        for (i = 0; i < rows; i++)
            for (j = 0; j < cols; j++)
            {
               neighbors = number_of_neighbors(i, j);
               switch(neighbors) 
               {
                   case 2:
                       new_grid[i * cols + j] = grid[i * cols + j];
                       break;
                   case 3:
                       new_grid[i * cols + j] = 1;
                       break;
                   case 0:
                   case 1:
                   default:
                        new_grid[i * cols + j] = 0;
                        break;
               }
            }

        // Make new grid to current grid for the next generation
        swap_ptr(&grid, &new_grid);
    }
#ifdef DEBUG
    printf("\nFinal Generation %d \n", MAXITER);
    print_grid();
#endif
}

void *bacteria_thread(void *arg)
{
    long id = (long)arg;
    int i, j, time, neighbors;
    
    int start_row = id * rows / NTHREADS;
    int end_row = (id + 1) * rows / NTHREADS;

    for (time = 0; time < MAXITER; time++)
    {
        for(i = start_row; i < end_row; i++)
            for(j = 0; j < cols; j++)
            {
                neighbors = number_of_neighbors(i, j);
                switch(neighbors) 
               {
                    case 2:
                        new_grid[i * cols + j] = grid[i * cols + j];
                        break;
                    case 3:
                        new_grid[i * cols + j] = 1;
                        break;
                    case 0:
                    case 1:
                    default:
                        new_grid[i * cols + j] = 0;
                        break;
               }
            }
        pthread_barrier_wait(&barrier);
        pthread_barrier_wait(&barrier2);
    }
    return NULL;
}
void parallel_bacteria(void)
{
    pthread_barrier_init(&barrier, NULL, NTHREADS+1);
    pthread_barrier_init(&barrier2, NULL, NTHREADS+1);
    
    for (int i = 0; i < NTHREADS; i++)
    {
        pthread_create(&threads[i], NULL, bacteria_thread, (void *)(long)i);
    }
    for(int i = 0; i < MAXITER; i++)
    {
#ifdef DEBUG
        printf("\nGeneration %d (Parallel)\n", i);
        print_grid();
#endif
        pthread_barrier_wait(&barrier);
        swap_ptr(&grid, &new_grid);
        pthread_barrier_wait(&barrier2);
    }
#ifdef DEBUG
    printf("\nFinal Generation %d (Parallel)\n", MAXITER);
    print_grid();
#endif
    for (int i = 0; i < NTHREADS; i++)
    {
        pthread_join(threads[i], NULL);
    }
    pthread_barrier_destroy(&barrier);
    pthread_barrier_destroy(&barrier2);
}

 
int main(int argc, char *argv[])
{
    struct timespec start, end;
    double serial, parallel;
    char output_filename[256];

    // Validate command line arguments
    if (argc < 4)
    {
        printf("Usage: %s <input_file> <num_generations> <num_threads>\n", argv[0]);
        printf("Example: %s bacteria1000.txt 250 10\n", argv[0]);
        return 1;
    }

    // Parse command line arguments
    MAXITER = atoi(argv[2]);
    NTHREADS = atoi(argv[3]);

    if (MAXITER <= 0 || NTHREADS <= 0)
    {
        printf("Error: num_generations and num_threads must be positive integers\n");
        return 1;
    }

    // Allocate threads array dynamically
    threads = (pthread_t *)malloc(NTHREADS * sizeof(pthread_t));
    if (!threads)
    {
        printf("Memory allocation error for threads\n");
        return 1;
    }


    read_from_file(argv[1]);
    printf("Initialize grid size Rows=%d, Cols=%d\n", rows, cols);

    printf("Start Serial with MAXITER=%d\n", MAXITER);
    clock_gettime(CLOCK_MONOTONIC, &start);
    serial_bacteria();
    clock_gettime(CLOCK_MONOTONIC, &end);
    serial = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Serial Time %lf \n", serial);
    
    strcpy(output_filename, argv[1]);
    *(strchr(output_filename, '.')) = '\0';
    strcat(output_filename, "_serial_out.txt");
    write_grid(output_filename);
    save_groundtruth(); // keep values from serial result as ground truth for later comparison

    printf("Initialize grid size Rows=%d, Cols=%d\n", rows, cols);
    read_from_file(argv[1]); // init again the same grid for parallel version

    printf("Start Paralel with NTHREADS=%d\n", NTHREADS);
    clock_gettime(CLOCK_MONOTONIC, &start);
    parallel_bacteria();
    clock_gettime(CLOCK_MONOTONIC, &end);
    parallel = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Parallel Time %lf  Speedup %lf \n", parallel, serial / parallel);
    
    strcpy(output_filename, argv[1]);
    *(strchr(output_filename, '.')) = '\0';
    strcat(output_filename, "_parallel_out.txt");
    write_grid(output_filename);
    
    if (!equal_groundtruth())
        printf("!!! Parallel version produces a different result! \n");
    else
        printf("Parallel version produced the same result \n");

    // Cleanup allocated memory
    free(threads);
    free(grid);
    free(new_grid);
    free(groundtruth);
    

    return 0;
}