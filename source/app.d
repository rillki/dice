module app;

void main(string[] args)
{
    import std;
    import dcv;

    auto image = imread("imgs/icon.png");
    imshow(image, "test");
}


