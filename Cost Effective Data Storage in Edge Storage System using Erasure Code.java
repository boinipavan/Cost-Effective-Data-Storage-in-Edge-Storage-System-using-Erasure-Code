import java.util.Scanner;

public class MatrixCreation {
    static int availrows;
    public static void main(String[] args) 
    {
        
        Scanner scanner = new Scanner(System.in);

        // Step 1: Take string as input
        System.out.print("Enter a string: ");
        String inputString = scanner.nextLine();

        // Step 2: Find ceil value of square root of length of string
        double ceilSqrt = Math.ceil(Math.sqrt(inputString.length()));

        // Step 3: Create n*n matrix (double data type) where n=ceil value of square root of length of string
        int n = (int) ceilSqrt;
        double[][] matrix = new double[n][n];

        // Step 4: Assign relevant ascii values of string characters sequentially
        int index = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (index < inputString.length()) {
                    matrix[i][j] = inputString.charAt(index);
                    index++;
                } else {
                    break;
                }
            }
            if (index >= inputString.length()) {
                break;
            }
        }

        // Step 5: If some part of matrix is unassigned then consider 'a' as parity fill remaining unassigned places with relevant ascii values
        char parityChar = 'a';
        while (index < n * n) {
            int row = index / n;
            int col = index % n;
            matrix[row][col] = parityChar;
            index++;
        
        }

        // Output the matrix
        System.out.println("Generated Matrix:");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
        //calling generateerasurematrix 
        double[][] erasurematrix=generateErasureMatrix();
         System.out.println("Generated erasureMatrix:");
        for (int i = 0; i < erasurematrix.length; i++) {
            for (int j = 0; j <erasurematrix[0].length; j++) {
                System.out.print(erasurematrix[i][j] + " ");
            }
            System.out.println();
        }
    
        double[][] codedmatrix = multiplyMatrices(erasurematrix, matrix);

        // Print the codedmatrixresult
        System.out.println("Resultant coded Matrix:");
        for (int i = 0; i < codedmatrix.length; i++) {
           for (int j = 0; j < codedmatrix[0].length; j++) {
               System.out.print(codedmatrix[i][j] + " ");
           }
            System.out.println();
        }
        //take input to know no of available rows
        System.out.print("enter avail rows:");
         availrows=scanner.nextInt();
         double[][] availableErasureMatrix=new double[availrows][];
        double[][] availableCodedMatrix = retrieveAvailableRows(codedmatrix,erasurematrix,availableErasureMatrix);

        // Output the available coded matrix
        System.out.println("Available Rows coded Matrix:");
        for (double[] row : availableCodedMatrix) {
            for (double element : row) {
                System.out.print(element + " ");
            }
            System.out.println();
        }
         // Output the available erasure matrix
        System.out.println("Available Rows erasure Matrix:");
        for (double[] row : availableErasureMatrix) {
            for (double element : row) {
                System.out.print(element + " ");
            }
            System.out.println();
        }
        //last matrix inverse operations
         // Find determinant
        double det = determinant(availableErasureMatrix);
        System.out.println("Determinant: " + det);

        // Find adjoint
        double[][] adj = adjoint(availableErasureMatrix);

        // Output adjoint
        System.out.println("Adjoint:");
        for (int i = 0; i < adj.length; i++) {
            for (int j = 0; j < adj[i].length; j++) {
                System.out.print(adj[i][j] + " ");
            }
            System.out.println();
        }
         double[][] result = multiplyMatrices(adj, availableCodedMatrix);

        // Output result
        System.out.println("Matrix Multiplication Result or actual matrix:");
        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[0].length; j++) {
                System.out.print(result[i][j]/det + " ");
            }
            System.out.println();
        }
           // Convert the integers to relevant strings based on ASCII values
        System.out.println("Strings based on ASCII values:");
        int k=0;
        for (int i = 0; i < result.length && k<inputString.length(); i++) {
            for (int j = 0; j <result[0].length && k<inputString.length(); j++) {
                System.out.print((char) result[i][j]/det + " ");
                k++;
            }
            System.out.println(); // Move to the next line after each row
        }
    }//static block 
    
       
    
      public static double[][] retrieveAvailableRows(double[][] originalMatrix,double[][] erasureMatrix,double[][] availableErasureMatrix) {
        Scanner scanner = new Scanner(System.in);
        
        //int numRows = originalMatrix.length;
        double[][] availableMatrix = new double[availrows][];
        int l=0;
        for (int i = 0; i <originalMatrix.length; i++) {
            System.out.print("Is row " + i + " available? (1 for available, 0 for not available): ");
            int choice = scanner.nextInt();
            if (choice == 1) {
                availableMatrix[l] = originalMatrix[i];
                availableErasureMatrix[l]=erasureMatrix[i];
                  l++;
            } 
          
        }

        return availableMatrix;
    }
    
    
     public static double[][] multiplyMatrices(double[][] matrix1, double[][] matrix2) {
        int rows1 = matrix1.length;
        int cols1 = matrix1[0].length;
        int rows2 = matrix2.length;
        int cols2 = matrix2[0].length;

        if (cols1 != rows2) {
            throw new IllegalArgumentException("Cannot multiply matrices of incompatible dimensions");
        }

        double[][] result = new double[rows1][cols2];

        for (int i = 0; i < rows1; i++) {
            for (int j = 0; j < cols2; j++) {
                for (int k = 0; k < cols1; k++) {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }

        return result;
    }
    //function to generate erasure matrix
    public static double[][] generateErasureMatrix()
    {
         Scanner scanner = new Scanner(System.in);

        // Input values for m and n
        System.out.print("Enter value for m (number of identity matrix rows): ");
        int m = scanner.nextInt();
        System.out.print("Enter value for n (number of parity matrix rows): ");
        int n = scanner.nextInt();

        // Create the matrix of size m + n rows
        double[][] matrix = new double[m + n][m];
        
        // Assign the first m rows as the identity matrix values
        for (int i = 0; i < m; i++) {
            matrix[i][i] = 1 ;
        }

        // Assign the remaining n rows as parity values starting from 1
        int parityValue = 1;
        for (int i = m; i < m + n; i++) {
            for (int j = 0; j < m; j++) {
                matrix[i][j] = parityValue;
                parityValue++;
            }
            
        }
        return matrix;

    }
    //matrix inverse
       // Function to calculate the cofactor of matrix[i][j]
    public static double[][] getCofactor(double[][] matrix, int i, int j) {
        int n = matrix.length;
        double[][] cofactor = new double[n - 1][n - 1];
        int row = 0, col = 0;
        for (int k = 0; k < n; k++) {
            for (int l = 0; l < n; l++) {
                if (k != i && l != j) {
                    cofactor[row][col++] = matrix[k][l];
                    if (col == n - 1) {
                        col = 0;
                        row++;
                    }
                }
            }
        }
        return cofactor;
    }

    // Function to find the determinant of a matrix
    public static double determinant(double[][] matrix) {
        int n = matrix.length;
        if (n == 1) {
            return matrix[0][0];
        }
        double det = 0;
        int sign = 1;
        for (int i = 0; i < n; i++) {
            double[][] cofactor = getCofactor(matrix, 0, i);
            det += sign * matrix[0][i] * determinant(cofactor);
            sign = -sign;
        }
        return det;
    }

    // Function to find the adjoint of a matrix
    public static double[][] adjoint(double[][] matrix) {
        int n = matrix.length;
        double[][] adjoint = new double[n][n];
        int sign;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                sign = ((i + j) % 2 == 0) ? 1 : -1;
                double[][] cofactor = getCofactor(matrix, i, j);
                adjoint[j][i] = sign * determinant(cofactor);
            }
        }
        return adjoint;
    }
        
  
}
