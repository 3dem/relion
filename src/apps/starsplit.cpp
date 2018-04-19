
#include <unistd.h>
#include <fstream>

#include <src/metadata_table.h>

void usage()
{
    std::cout << "usage: relion_starsplit <.star> <label> <cat1> <cat2> ...\n";
}

int main(int argc, char *argv[])
{

    if (argc < 5)
    {
        usage();
        return 666;
    }

    std::string fn = argv[1];

    MetaDataTable mdt;
    mdt.read(fn);

    std::string labelName = argv[2];
    if (labelName[0] == '_')
    {
        labelName = labelName.substr(1);
    }

    std::cout << "label = " << labelName << "\n";

    EMDLabel label = EMDL::str2Label(labelName);

    if (label == EMDL_UNDEFINED)
    {
        std::cout << "invalid label: " << labelName << "\n";
        return 666;
    }

    int catNum = argc - 3;

    std::vector<std::string> cats(catNum);

    for (int i = 0; i < catNum; i++)
    {
        cats[i] = argv[i+3];
    }

    std::vector<MetaDataTable> tables(catNum);

    for (int i = 0; i < mdt.numberOfObjects(); i++)
    {
        std::string val;
        mdt.getValueToString(label, val, i);

        for (int j = 0; j < catNum; j++)
        {
            if (val.find(cats[j]) != std::string::npos)
            {
                //std::cout << "found " << cats[j] << " in " << val << " at " << val.find(cats[j]) << "\n";
                //std::cout << val.find(cats[j]) << " != " << std::string::npos << "\n";

                tables[j].addObject(mdt.getObject(i));
                break;
            }
        }
    }

    std::string base = fn.substr(0, fn.length()-5);

    for (int i = 0; i < catNum; i++)
    {
        tables[i].write(base+"_"+cats[i]+".star");
    }

    return 0;
}
