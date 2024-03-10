import java.util.Scanner;
import javax.crypto.Cipher;
import javax.crypto.KeyGenerator;
import javax.crypto.SecretKey;
import javax.crypto.spec.SecretKeySpec;
import java.util.Arrays;
import java.util.Base64;

public class MatrixCreation {
    static int availrows,n;
    static double ceilSqrt;
    public static void main(String[] args) throws Exception
    {
        
        Scanner scanner = new Scanner(System.in);

        // Step 1: Take string as input
        System.out.print("Enter a string(Data): ");
        String inputString = scanner.nextLine();

        // Step 2: Find ceil value of square root of length of string
        ceilSqrt = Math.ceil(Math.sqrt(inputString.length()));

        // Step 3: Create n*n matrix (double data type) where n=ceil value of square root of length of string
         n = (int)ceilSqrt;
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
        System.out.println("Generated matrix:");
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
        //perform aes des encryption
        SecretKeySpec aesKeySpec = generateAESKey();
        SecretKeySpec desKeySpec = generateDESKey();
        String[] encryptedRows = new String[codedmatrix.length];
        for(int i=0;i<codedmatrix.length;i++)
        {
                  StringBuilder sb1 = new StringBuilder();
                  for (double row : codedmatrix[i]) {
                      
                          char c = (char) row;
                          sb1.append(c);
                      }
                  String s=sb1.toString();
                  String encryptedWithAES = encryptWithAES(s, aesKeySpec);
                  encryptedRows[i] = encryptWithDES(encryptedWithAES, desKeySpec);
        }
        
        
        //take input to know no of available rows
        	System.out.print("enter number of available servers(rows):");
        	availrows=scanner.nextInt();
         if(availrows<n)
         {
        	 System.out.print("Number of servers failed is greater than parity servers.Hence,we cant recover data.sorry..");
        	 System.exit(0);
         }
         double[][] availableErasureMatrix=new double[matrix.length][];
     	
        
         String[] availableCodedMatrix = retrieveAvailableRows(encryptedRows,erasurematrix,availableErasureMatrix);//encrypted
        
        
        //perform des aes decryption
        
            double[][] availableEncodedMatrix=new double[matrix.length][matrix.length];
        	for(int x=0;x<availableCodedMatrix.length;x++)
        	{
        		String s=availableCodedMatrix[x];
        		String decryptedWithDES = decryptWithDES(s, desKeySpec);
        		String str = decryptWithAES(decryptedWithDES, aesKeySpec);
            for (int i = 0; i < str.length(); i++) {
                char c = str.charAt(i);
             
                availableEncodedMatrix[x][i] = (double) (int) c; // Store the ASCII value as double
            }
        	}

        // Output the available coded matrix
        System.out.println("Available Rows coded Matrix:");
        for (double[] row : availableEncodedMatrix) {
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
        if(det==0)
        {
        System.out.println("Determinant is zero i.e;matrix is scalar please consider differrent parity values in erasure matrix");
        }

        // Find adjoint
        double[][] adj = adjoint(availableErasureMatrix);

        // Output adjoint
        /*
        System.out.println("Adjoint:");
        for (int i = 0; i < adj.length; i++) {
            for (int j = 0; j < adj[i].length; j++) {
                System.out.print(adj[i][j] + " ");
            }
            System.out.println();
        }*/
         double[][] result = multiplyMatrices(adj, availableEncodedMatrix);

        // Output result
        System.out.println("Matrix Multiplication Result or actual matrix after data re-construction");
        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[0].length; j++) {
                System.out.print(result[i][j]/det + " ");
            }
            System.out.println();
        }
           // Convert the integers to relevant strings based on ASCII values
            for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[0].length; j++) {
                result[i][j]=result[i][j]/det;
            }
            }
        System.out.print("Data Stored is:");
        int k=0;
        for (int i = 0; i < result.length && k<inputString.length(); i++) {
            for (int j = 0; j <result[0].length && k<inputString.length(); j++) {
                System.out.print((char) result[i][j]);
                k++;
            }
            
        }
    }//static block 
    
       
    
      public static String[] retrieveAvailableRows(String[] originalMatrix,double[][] erasureMatrix,double[][] availableErasureMatrix) {
        Scanner scanner = new Scanner(System.in);
        
        //int numRows = originalMatrix.length;
        String[] availableMatrix = new String[n];
        int l=0;
        for (int i = 0; i <originalMatrix.length && l<availableErasureMatrix.length; i++) {
            System.out.print("Is server(row) " + i + " available? (1 for available, 0 for not available): ");
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
        //System.out.print("Enter value for m (number of identity matrix rows): ");
        int m = (int)ceilSqrt;
        System.out.print("Enter number of parity servers(rows): ");
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
    public static SecretKeySpec generateAESKey() throws Exception {
        KeyGenerator keyGen = KeyGenerator.getInstance("AES");
        keyGen.init(128); // AES key length can be 128, 192, or 256 bits
        SecretKey secretKey = keyGen.generateKey();
        return new SecretKeySpec(secretKey.getEncoded(), "AES");
    }

    public static SecretKeySpec generateDESKey() throws Exception {
        KeyGenerator keyGen = KeyGenerator.getInstance("DES");
        keyGen.init(56); // DES key length is 56 bits
        SecretKey secretKey = keyGen.generateKey();
        return new SecretKeySpec(secretKey.getEncoded(), "DES");
    }

    public static String encryptWithAES(String message, SecretKeySpec key) throws Exception {
        Cipher aesCipher = Cipher.getInstance("AES");
        aesCipher.init(Cipher.ENCRYPT_MODE, key);
        byte[] encryptedBytes = aesCipher.doFinal(message.getBytes());
        return Base64.getEncoder().encodeToString(encryptedBytes);
    }

    public static String encryptWithDES(String message, SecretKeySpec key) throws Exception {
        Cipher desCipher = Cipher.getInstance("DES");
        desCipher.init(Cipher.ENCRYPT_MODE, key);
        byte[] encryptedBytes = desCipher.doFinal(message.getBytes());
        return Base64.getEncoder().encodeToString(encryptedBytes);
    }

    public static String decryptWithDES(String encryptedText, SecretKeySpec key) throws Exception {
        byte[] encryptedBytes = Base64.getDecoder().decode(encryptedText);
        Cipher desCipher = Cipher.getInstance("DES");
        desCipher.init(Cipher.DECRYPT_MODE, key);
        byte[] decryptedBytes = desCipher.doFinal(encryptedBytes);
        return new String(decryptedBytes);
    }

    public static String decryptWithAES(String encryptedText, SecretKeySpec key) throws Exception {
        byte[] encryptedBytes = Base64.getDecoder().decode(encryptedText);
        Cipher aesCipher = Cipher.getInstance("AES");
        aesCipher.init(Cipher.DECRYPT_MODE, key);
        byte[] decryptedBytes = aesCipher.doFinal(encryptedBytes);
        return new String(decryptedBytes);
    }
        
  
}
