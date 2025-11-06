function x = Gauss_Seidel(A,B)

tolerance = 1e-4;
x = zeros(size(A,1),1);
error = 10;

while error>tolerance
    x_old = x;
    for i = 1:size(A,1)
        sum = 0;
        for j =  1:i-1
            sum = sum + A(i,j)*x(j);
        end  

        for j = i+1:size(A,1)
            sum = sum + A(i,j)*x_old(j);
        end  
        x(i) = (1/(A(i,i)))*(B(i)-sum);

    end
    error = norm(x - x_old,2);
end
disp(x)

end