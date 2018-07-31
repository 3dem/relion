#include <memory>
#include <iostream>

class Wrapper
{
    public:

        Wrapper()
        {
            payload = std::shared_ptr<Payload>(new Payload);
        }


        class Payload
        {
            public:

                Payload(){}

                ~Payload()
                {
                    std::cout << "deleting payload\n";
                }
        };

        std::shared_ptr<Payload> payload;
};


int main(int argc, char *argv[])
{
    Wrapper v0;

    std::cout << "0\n";
    {
        Wrapper w0;

        std::cout << "1\n";

        {
            Wrapper w1 = w0;
            std::cout << "2\n";
        }

        std::cout << "3\n";
    }

    std::cout << "4\n";
}
