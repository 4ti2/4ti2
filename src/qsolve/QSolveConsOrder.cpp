#include "qsolve/QSolveConsOrder.h"

using namespace _4ti2_;

ConsOrder::ConsOrder()
{
    // The default constraint order.
    set_constraint_order(MAXCUTOFF);
}

ConsOrder::ConsOrder(QSolveConsOrder o)
{
    set_constraint_order(o);
}

void
ConsOrder::set_constraint_order(QSolveConsOrder o)
{
    order = o;
    switch (o)
    {
    case MININDEX:
        compare = &minindex;
        break;
    case MAXCUTOFF:
        compare = &maxcutoff;
        break;
    case MINCUTOFF:
        compare = &mincutoff;
        break;
    case MAXINTER:
        compare = &maxinter;
        break;
    default:
        std::cerr << "Unrecognized constraint order.\n";
        exit(1);
        break;
    }
    circuit_compare = &maxinter;
}

QSolveConsOrder
ConsOrder::get_constraint_order() const
{
    return order;
}

bool
ConsOrder::minindex(Index next_pos_count, Index next_neg_count, Index next_zero_count,
                Index pos_count, Index neg_count, Index zero_count)
{
    return false;
}

bool
ConsOrder::maxinter(Index next_pos_count, Index next_neg_count, Index next_zero_count,
                Index pos_count, Index neg_count, Index zero_count)
{
    return (zero_count > next_zero_count);
}

bool
ConsOrder::maxcutoff(Index next_pos_count, Index next_neg_count, Index next_zero_count,
                Index pos_count, Index neg_count, Index zero_count)
{
    return (neg_count > next_neg_count);
}

bool
ConsOrder::mincutoff(Index next_pos_count, Index next_neg_count, Index next_zero_count,
                Index pos_count, Index neg_count, Index zero_count)
{
    return (neg_count < next_neg_count);
}

