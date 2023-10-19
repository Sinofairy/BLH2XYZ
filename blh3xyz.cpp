#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <eigen3/Eigen/Eigen>
using namespace std;

const double epsilon = 0.000000000000001;
const double pi = 3.14159265358979323846;
const double d2r = pi / 180;
const double r2d = 180 / pi;

const double a = 6378137.0;                //椭球长半轴
const double f_inverse = 298.257223563;                        //扁率倒数
const double b = a - a / f_inverse;
//const double b = 6356752.314245;                        //椭球短半轴

const double e = sqrt(a * a - b * b) / a;


void Blh2Xyz(double &x, double &y, double &z)  //经纬高转地心地固坐标系
{
        double L = x * d2r;
        double B = y * d2r;
        double H = z;

        double N = a / sqrt(1 - e * e * sin(B) * sin(B));
        x = (N + H) * cos(B) * cos(L);
        y = (N + H) * cos(B) * sin(L);
        z = (N * (1 - e * e) + H) * sin(B);
}


void CalEcef2Enu(Eigen::Vector3d& topocentricOrigin, Eigen::Matrix4d& resultMat)
{
        double rzAngle = -(topocentricOrigin.x() * d2r + pi / 2);
        Eigen::AngleAxisd rzAngleAxis(rzAngle, Eigen::Vector3d(0, 0, 1));
        Eigen::Matrix3d rZ = rzAngleAxis.matrix();

        double rxAngle = -(pi / 2 - topocentricOrigin.y() * d2r);
        Eigen::AngleAxisd rxAngleAxis(rxAngle, Eigen::Vector3d(1, 0, 0));
        Eigen::Matrix3d rX = rxAngleAxis.matrix();

        Eigen::Matrix4d rotation;
        rotation.setIdentity();
        rotation.block<3, 3>(0, 0) = (rX * rZ);
        //cout << rotation << endl;
                                
        double tx = topocentricOrigin.x();
        double ty = topocentricOrigin.y();
        double tz = topocentricOrigin.z();
        Blh2Xyz(tx, ty, tz);
        Eigen::Matrix4d translation;
        translation.setIdentity();
        translation(0, 3) = -tx;
        translation(1, 3) = -ty;
        translation(2, 3) = -tz;
        
        resultMat = rotation * translation;
}

int main() {

    vector<vector<string>> data; // 定义一个二维 vector 用来存储数据

    ifstream in_file("gps.txt"); // 打开文件
    if (!in_file) {
        cerr << "无法打开文件！" << endl;
        return 1;
    }

    string line;
    while (getline(in_file, line)) { // 逐行读取文件内容
        vector<string> row; // 定义一个 vector 用来存储每一行的数据

        stringstream ss(line);
        string field;
        while (getline(ss, field, ' ')) { // 按照空格分割每一行的内容
            row.push_back(field); // 将分割后的内容存储到 vector 中
        }

        data.push_back(row); // 将每一行的数据存储到二维 vector 中
    }

    in_file.close(); // 关闭文件

    // 输出二维 vector 中的数据
    // for (const auto& row : data) {
    //     for (const auto& field : row) {
    //         cout << field << " ";
    //     }
    //     cout << endl;
    // }
    // double TimeStamp = 0.0;
    double L = 0.0;
    double B = 0.0;
    double H = 0.0;
    double Yaw = 0.0;
    
    double TimeStamp = std::stof(data[0][0]);
    double L0 = std::stod(data[0][1]);
    double B0 = std::stod(data[0][2]);
    double H0 = std::stod(data[0][3]);
    double Yaw0 = std::stod(data[0][4]);
    double Delta_Yaw = 0.0;
    Eigen::Vector3d topocentricOrigin(L0, B0, H0);
    Eigen::Matrix4d wolrd2localMatrix;
    CalEcef2Enu(topocentricOrigin, wolrd2localMatrix);        
    cout << "地心转站心矩阵：" << endl;
    cout << wolrd2localMatrix << endl<<endl;
    Blh2Xyz(L0, B0, H0);
    // cout << "ECEF坐标（世界坐标）：";
    Eigen::Vector4d xyz(L0, B0, H0, 1);
    // cout << xyz << endl;

    // cout << "ENU坐标（局部坐标）原点：";
    Eigen::Vector4d enu = wolrd2localMatrix * xyz;
    cout << TimeStamp << " " << enu[0] << " " << enu[1] << " " << enu[2] << " "<< Delta_Yaw << endl;  

    for(int i = 1; i < data.size(); i++)
    {
        TimeStamp = std::stof(data[i][0]);
        L = std::stod(data[i][1]);
        B = std::stod(data[i][2]);
        H = std::stod(data[i][3]);
        Yaw = std::stod(data[i][4]);
        if(Yaw - 0.0 < -epsilon)
        {
            Yaw += 360;
        }

        Delta_Yaw = Yaw - Yaw0;
        
        if(Delta_Yaw > 180)
        {
            Delta_Yaw = -(360 - Delta_Yaw);
        }
        if(Delta_Yaw != 0.0){
            Delta_Yaw = -Delta_Yaw;
        }
        
        // Eigen::Vector3d topocentricOrigin(L, B, H);
        // Eigen::Matrix4d wolrd2localMatrix;
        // CalEcef2Enu(topocentricOrigin, wolrd2localMatrix);        
        // cout << "地心转站心矩阵：" << endl;
        // cout << wolrd2localMatrix << endl<<endl;

        Blh2Xyz(L, B, H);

        // cout << "ECEF坐标（世界坐标）：";
        Eigen::Vector4d xyz(L, B, H, 1);
        // cout << xyz << endl;

        // cout << "ENU坐标（局部坐标）：";
        Eigen::Vector4d enu = wolrd2localMatrix * xyz;
        //cout << TimeStamp << " " << std::setw(10) << enu[0] << " " << enu[1] << " " << enu[2] << " "<< Delta_Yaw << endl;  
        cout << TimeStamp << " " << enu[0] << " " << enu[1] << " " << enu[2] << " "<< Delta_Yaw << endl;  
    }


    return 0;
}