//
// Created by Ronald on 20/02/2020.
//

template <typename d_type, typename vec_type>
class respiration_co2{
public:

    respiration_co2(pear::component<d_type, vec_type> & co2, pear::component<d_type, vec_type> & o2, pear::grid<d_type> & grid,
                    d_type p1,
                    d_type p2,
                    d_type p3
    )
            : o2_(o2)
            , co2_(co2)
            , grid_(grid)
            , p1_(p1)
            , p2_(p2)
            , p3_(p3)
    {
        std::cout<<"Component "<<co2_.name()<<" reacts with "<<o2_.name()<<std::endl;
    }

private:
    pear::component<d_type, vec_type> & o2_;
    pear::component<d_type, vec_type> & co2_;
    pear::grid<d_type> & grid_;
    d_type p1_, p2_, p3_;

};
