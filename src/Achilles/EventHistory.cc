#include "Achilles/EventHistory.hh"

using achilles::EventHistoryNode;
using achilles::EventHistory;

void EventHistory::AddParticleIn(size_t idx, const Particle &part) {
    m_history[idx] -> AddIncoming(part);
}

void EventHistory::AddParticleOut(size_t idx, const Particle &part) {
    m_history[idx] -> AddOutgoing(part);
}

std::vector<EventHistoryNode*> EventHistory::Children(size_t idx) const {
    std::vector<EventHistoryNode*> children;
    for(const auto &part : m_history[idx] -> ParticlesOut()) {
        auto node = FindNodeIn(part);
        if(node) children.push_back(node);
    }
    return children;
}

std::vector<EventHistoryNode*> EventHistory::Parents(size_t idx) const {
    std::vector<EventHistoryNode*> parents;
    for(const auto &part : m_history[idx] -> ParticlesIn()) {
        auto node = FindNodeOut(part);
        if(node) parents.push_back(node);
    }
    return parents;
}

std::vector<EventHistoryNode*> EventHistory::Children(EventHistoryNode *node) const {
    auto it = std::find_if(m_history.cbegin(), m_history.cend(), [&](const std::unique_ptr<EventHistoryNode> &p) {
        return p.get() == node;
    });
    if(*it == nullptr) return {};
    return Children((*it)->Index());
}

std::vector<EventHistoryNode*> EventHistory::Parents(EventHistoryNode *node) const {
    auto it = std::find_if(m_history.cbegin(), m_history.cend(), [&](const std::unique_ptr<EventHistoryNode> &p) {
        return p.get() == node;
    });
    if(it == m_history.end()) return {};
    return Parents((*it)->Index());
}

EventHistoryNode* EventHistory::Node(size_t idx) const {
    for(const auto &node : m_history) {
        if(node -> Index() == idx) {
            return node.get();
        }
    }
    return nullptr;
}

EventHistoryNode* EventHistory::GetUniqueNode(StatusCode status) const {
    size_t n_nodes{};
    EventHistoryNode *result = nullptr;
    for(const auto &node : m_history) {
        if(node -> Status() == status) {
            result = node.get();
            n_nodes++;
        }
        if(n_nodes > 1) {
            throw std::runtime_error(fmt::format("EventHistory: Only one {} node is allowed!",
                                     ToString(status)));
        }
    }
    return result;
}

EventHistoryNode* EventHistory::FindNode(bool incoming, const Particle &part) const {
    for(const auto &node : m_history) {
        std::vector<Particle> particles;

        if(incoming) {
            if(node -> HasIncoming(part)) return node.get();
        } else {
            if(node -> HasOutgoing(part)) return node.get();
        }
    }
    return nullptr;
}
