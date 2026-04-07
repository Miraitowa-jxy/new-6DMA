function Q = rotation_matrix(psi)
Qz =[cos(psi(1)), -sin(psi(1)), 0;
    sin(psi(1)), cos(psi(1)), 0;
    0, 0, 1] ; %绕z轴的旋转矩阵
Qy =[cos(psi(2)), 0, sin(psi(2));
    0, 1, 0;
    -sin(psi(2)), 0, cos(psi(2))]; 
Qx = [1, 0, 0
    0, cos(psi(3)), -sin(psi(3));
    0, sin(psi(3)), cos(psi(3))];
Q = Qz * Qy * Qx;
end

