//
// Created by dell on 2022/7/13.
//

#ifndef TEST_KALMAN_H
#define TEST_KALMAN_H

typedef float              float32_t;
Matrix m_StatePost;
void init(int32_t stateNum, int32_t measurementNum)
{
    m_StatePost = Matrix(stateNum ,1);
    m_StatePost.zero();
}


void setState(const SlotTrackBox& state)
{
    m_StatePost.val[0][0]  = 1/*state.pos1.lateral*/;            m_StatePost.val[2][0]  = 1;
    m_StatePost.val[1][0]  = 5/*state.pos1.longitudinal*/;       m_StatePost.val[3][0]  = 1;

    m_StatePost.val[4][0]  = 2/*state.pos2.lateral*/;            m_StatePost.val[6][0]  = 1;
    m_StatePost.val[5][0]  = 5/*state.pos2.longitudinal*/;       m_StatePost.val[7][0]  = 1;

    m_StatePost.val[8][0]  = 2/*state.pos3.lateral*/;            m_StatePost.val[10][0] = 1;
    m_StatePost.val[9][0]  = 1/*state.pos3.longitudinal*/;       m_StatePost.val[11][0] = 1;

    m_StatePost.val[12][0] = 1/*state.pos4.lateral*/;            m_StatePost.val[14][0] = 1;
    m_StatePost.val[13][0] = 1/*state.pos4.longitudinal*/;       m_StatePost.val[15][0] = 1;

    float32_t sin       =   sinf(yawRate);
    float32_t cos       =   cosf(yawRate);
    float32_t deltaX    =   distance * sinf(yawRate);
    float32_t deltaY    =   distance * cosf(yawRate);
    Matrix    angle     =   Matrix(16, 16);

    angle.zero();
    angle.val[0][0]     =   cos;        angle.val[0][1]    =  -sin;
    angle.val[1][0]     =   sin;        angle.val[1][1]    =   cos;
    angle.val[2][2]     =   1;          angle.val[3][3]    =     1;
    angle.val[4][4]     =   cos;        angle.val[4][5]    =  -sin;
    angle.val[5][4]     =   sin;        angle.val[5][5]    =   cos;
    angle.val[6][6]     =   1;          angle.val[7][7]    =     1;
    angle.val[8][8]     =   cos;        angle.val[8][9]    =  -sin;
    angle.val[9][8]     =   sin;        angle.val[9][9]    =   cos;
    angle.val[10][10]   =   1;          angle.val[11][11]  =     1;
    angle.val[12][12]   =   cos;        angle.val[12][13]  =  -sin;
    angle.val[13][12]   =   sin;        angle.val[13][13]  =   cos;
    angle.val[14][14]   =   1;          angle.val[15][15]  =     1;

    m_F.zero();
    m_F.val[0][0]       =   1.0;        m_F.val[8][8]      =   1.0;
    m_F.val[1][1]       =   1.0;        m_F.val[9][9]      =   1.0;
    m_F.val[2][2]       =   1.0;        m_F.val[10][10]    =   1.0;
    m_F.val[3][3]       =   1.0;        m_F.val[11][11]    =   1.0;
    m_F.val[4][4]       =   1.0;        m_F.val[12][12]    =   1.0;
    m_F.val[5][5]       =   1.0;        m_F.val[13][13]    =   1.0;
    m_F.val[6][6]       =   1.0;        m_F.val[14][14]    =   1.0;
    m_F.val[7][7]       =   1.0;        m_F.val[15][15]    =   1.0;

    float32_t sinR      =   sinf(-1.0f * yawRate);
    float32_t cosR      =   cosf(-1.0f * yawRate);
    m_F.val[0][2]       =   deltaX * sinR;     m_F.val[8][10]     = deltaX * sinR;
    m_F.val[1][3]       =   deltaY * cosR;     m_F.val[9][11]     = deltaY * cosR;
    m_F.val[4][6]       =   deltaX * sinR;     m_F.val[12][14]    = deltaX * sinR;
    m_F.val[5][7]       =   deltaY * cosR;     m_F.val[13][15]    = deltaY * cosR;

    m_F = angle * m_F;

}
#endif //TEST_KALMAN_H
