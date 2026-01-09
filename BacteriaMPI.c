
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include <stdint.h>

// Global variables for command line arguments
int MAXITER;
int NPROCS;

//#define DEBUG

// MPI variables
int rank, size;

// Local grid for each process (includes ghost rows)
char *local_grid;    // size: (local_rows + 2) * cols, includes top and bottom ghost rows
char *local_new_grid;

// Full grid only for rank 0
char *grid;
char *new_grid;
char *groundtruth;
int rows, cols;

// Local rows per process
int local_rows;
int start_row, end_row;
 
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
    
    // Validate dimensions
    if (rows <= 0 || cols <= 0)
    {
        printf("Error: dimensions must be positive (rows=%d, cols=%d)\n", rows, cols);
        fclose(file);
        return 0;
    }
    
    // Check for integer overflow before multiplication
    if (cols > 0 && rows > SIZE_MAX / cols)
    {
        printf("Error: grid size too large (would cause overflow)\n");
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
    // Check for integer overflow (already validated in read_from_file, but being safe)
    if (cols > 0 && rows > SIZE_MAX / cols)
    {
        printf("Error: groundtruth size too large (would cause overflow)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    groundtruth = (char *)malloc(rows * cols * sizeof(char));
    if (!groundtruth)
    {
        printf("Memory allocation error for groundtruth result\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
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

int number_of_neighbors(char *g, int i, int j, int max_rows, int max_cols)
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
            
            if (ni >= 0 && ni < max_rows && nj >= 0 && nj < max_cols)
            {
                count += g[ni * max_cols + nj];
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
               neighbors = number_of_neighbors(grid, i, j, rows, cols);
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

void mpi_bacteria(void)
{
    int i, j, time, neighbors;
    
    // Calculate total local rows including ghost rows
    // Ghost row at top (index 0) and ghost row at bottom (index local_rows+1)
    int total_local_rows = local_rows + 2;
    
    // Check for integer overflow before allocation
    if (cols > 0 && total_local_rows > SIZE_MAX / cols)
    {
        printf("Rank %d: Local grid size too large (would cause overflow)\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // Allocate local grids with ghost rows
    local_grid = (char *)malloc(total_local_rows * cols * sizeof(char));
    local_new_grid = (char *)malloc(total_local_rows * cols * sizeof(char));
    
    if (!local_grid || !local_new_grid)
    {
        printf("Rank %d: Memory allocation error for local grids\n", rank);
        free(local_grid);
        free(local_new_grid);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // Initialize ghost rows to 0
    memset(local_grid, 0, total_local_rows * cols * sizeof(char));
    memset(local_new_grid, 0, total_local_rows * cols * sizeof(char));

    // Distribute rows from rank 0 to all processes
    int *sendcounts = NULL;
    int *displs = NULL;
    
    if (rank == 0)
    {
        sendcounts = (int *)malloc(size * sizeof(int));
        displs = (int *)malloc(size * sizeof(int));
        
        if (!sendcounts || !displs)
        {
            printf("Rank %d: Memory allocation error for sendcounts/displs\n", rank);
            free(sendcounts);
            free(displs);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        for (int p = 0; p < size; p++)
        {
            int p_start = p * rows / size;
            int p_end = (p + 1) * rows / size;
            sendcounts[p] = (p_end - p_start) * cols;
            displs[p] = p_start * cols;
        }
    }
    
    // Scatter the initial grid (into row 1 onwards, skipping ghost row at index 0)
    MPI_Scatterv(grid, sendcounts, displs, MPI_CHAR,
                 &local_grid[cols], local_rows * cols, MPI_CHAR,
                 0, MPI_COMM_WORLD);
    
    // Main evolution loop
    for (time = 0; time < MAXITER; time++)
    {   
        // Exchange border with upper neighbor (rank-1)
        if (rank > 0)
        {
            // Send my first real row (index 1) to rank-1
            MPI_Send(&local_grid[1 * cols], cols, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD);
            // Receive from rank-1 into my top ghost row (index 0)
            MPI_Recv(&local_grid[0 * cols], cols, MPI_CHAR, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        // Exchange border with lower neighbor (rank+1)
        if (rank < size - 1)
        {
            // Receive from rank+1 into my bottom ghost row (index local_rows+1)
            MPI_Recv(&local_grid[(local_rows + 1) * cols], cols, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // Send my last real row (index local_rows) to rank+1
            MPI_Send(&local_grid[local_rows * cols], cols, MPI_CHAR, rank + 1, 1, MPI_COMM_WORLD);
        }
        
        // Compute new generation for local rows (indices 1 to local_rows)
        for (i = 1; i <= local_rows; i++)
        {
            for (j = 0; j < cols; j++)
            {
                neighbors = number_of_neighbors(local_grid, i, j, total_local_rows, cols);
                switch(neighbors) 
                {
                    case 2:
                        local_new_grid[i * cols + j] = local_grid[i * cols + j];
                        break;
                    case 3:
                        local_new_grid[i * cols + j] = 1;
                        break;
                    case 0:
                    case 1:
                    default:
                        local_new_grid[i * cols + j] = 0;
                        break;
                }
            }
        }
        
        // Swap local grids
        swap_ptr(&local_grid, &local_new_grid);
    }
    
    // Gather results back to rank 0 (from row 1 onwards, skipping ghost row)
    MPI_Gatherv(&local_grid[cols], local_rows * cols, MPI_CHAR,
                grid, sendcounts, displs, MPI_CHAR,
                0, MPI_COMM_WORLD);
    
    // Cleanup
    if (rank == 0)
    {
        free(sendcounts);
        free(displs);
    }
    free(local_grid);
    free(local_new_grid);
}

 
int main(int argc, char *argv[])
{
    struct timespec start, end;
    double serial, parallel;
    char output_filename[256];

    // Initialize MPI first
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    NPROCS = size;

    // Validate command line arguments
    if (argc < 3)
    {
        if (rank == 0)
        {
            printf("Usage: mpirun -np <num_processes> %s <input_file> <num_generations>\n", argv[0]);
            printf("Example: mpirun -np 4 %s bacteria1000.txt 250\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    // Parse command line arguments
    MAXITER = atoi(argv[2]);

    if (MAXITER <= 0)
    {
        if (rank == 0)
            printf("Error: num_generations must be positive integer\n");
        MPI_Finalize();
        return 1;
    }
    
    // Check filename length (22 = length of "_parallel_out.txt" suffix)
    if (strlen(argv[1]) > 256 - 22)
    {
        if (rank == 0)
            printf("Error: filename too long (max %d characters)\n", 256 - 22);
        MPI_Finalize();
        return 1;
    }

    // Only rank 0 reads the file and runs serial version
    if (rank == 0)
    {
        if (!read_from_file(argv[1]))
        {
            printf("Failed to read input file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        printf("Initialize grid size Rows=%d, Cols=%d\n", rows, cols);

        printf("Start Serial with MAXITER=%d\n", MAXITER);
        clock_gettime(CLOCK_MONOTONIC, &start);
        serial_bacteria();
        clock_gettime(CLOCK_MONOTONIC, &end);
        serial = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        printf("Serial Time %lf \n", serial);
        
        strcpy(output_filename, argv[1]);
        char *dot = strchr(output_filename, '.');
        if (dot != NULL)
        {
            *dot = '\0';
        }
        strcat(output_filename, "_serial_out.txt");
        write_grid(output_filename);
        save_groundtruth(); // keep values from serial result as ground truth for later comparison

        printf("Initialize grid size Rows=%d, Cols=%d\n", rows, cols);
        free(grid);
        free(new_grid);
        if (!read_from_file(argv[1]))
        {
            printf("Failed to re-read input file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    if (rank == 0)
    {
        printf("Start Parallel with NPROCS=%d\n", NPROCS);
        clock_gettime(CLOCK_MONOTONIC, &start);
    }

    // Broadcast dimensions to all processes
    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Calculate local rows for each process
    //This is the formula that helps to distribute the rows equally between the processes
    //For example if we have 10 rows and 4 processes, the rows will be distributed like this:
    // 0: 0-1
    // 1: 2-4
    // 2: 5-6
    // 3: 7-9
    start_row = rank * rows / size;
    end_row = (rank + 1) * rows / size;
    local_rows = end_row - start_row;
    
    
    MPI_Barrier(MPI_COMM_WORLD);

    mpi_bacteria();
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Only rank 0 writes output and compares results
    if (rank == 0)
    {
        clock_gettime(CLOCK_MONOTONIC, &end);
        parallel = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        printf("Parallel Time %lf  Speedup %lf \n", parallel, serial / parallel);
        
        strcpy(output_filename, argv[1]);
        char *dot = strchr(output_filename, '.');
        if (dot != NULL)
        {
            *dot = '\0';
        }
        strcat(output_filename, "_parallel_out.txt");
        write_grid(output_filename);
        
        if (!equal_groundtruth())
            printf("!!! Parallel version produces a different result! \n");
        else
            printf("Parallel version produced the same result \n");

        // Cleanup allocated memory
        free(grid);
        free(new_grid);
        free(groundtruth);
    }
    
    MPI_Finalize();
    return 0;
}

/*
└─$ mpiexec -np 3 ./prog bacteria1000.txt 250
Initialize grid size Rows=1000, Cols=1000
Start Serial with MAXITER=250
Serial Time 8.104392
Grid saved to bacteria1000_serial_out.txt
Initialize grid size Rows=1000, Cols=1000
Start Parallel with NPROCS=3
Parallel Time 3.269405  Speedup 2.478858
Grid saved to bacteria1000_parallel_out.txt
Parallel version produced the same result
*/