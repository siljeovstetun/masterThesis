class People
{
  public:
    int condition; // 0 = healty, 1 = sick, 2 = imune
    int infectedTime;
    double xPosition;
    double yPosition;
    double angle;
    double v0;
    double dr;
    int numberOfInfected;

    People();

    People(int L, int i, int N, double v, double d);

    void setCondition(int condition, int time);

    void setPosition(double x, double y, double a);

    void setNumberOfInfected();

    void setv0(double v);

    void setDr(double d);

    int getCondition();

    double getxPosition();

    double getyPosition();

    double getAngle();

    double getv0();

    double getDr();

    int getTime();

    int getNumberOfInfected();
};
