function ROTL(a::T, b) where {T}
    return (a << b) | (a >> (32 - b))
end

function QR(block::Vector{T}, a, b, c, d) where {T}
    block[a] += block[b]
    block[d] ⊻= block[a]
    block[d] = ROTL(block[d], 16)

    block[c] += block[d]
    block[b] ⊻= block[c]
    block[b] = ROTL(block[b], 12)

    block[a] += block[b]
    block[d] ⊻= block[a]
    block[d] = ROTL(block[d], 8)

    block[c] += block[d]
    block[b] ⊻= block[c]
    block[b] = ROTL(block[b], 7)
end

function encrypt(key::Vector{T}, nonce::Vector{T}, counter::T) where {T}
    state = zeros(T, 16)
    state[1] = 0x61707865
    state[2] = 0x3320646e
    state[3] = 0x79622d32
    state[4] = 0x6b206574

    state[5:12] = copy(key) # Potentially do not need
    state[13] = counter
    state[14:16] = copy(nonce) # Potentially do not need, also things like views which may be helpful to look at
    print(state)
    initialState = copy(state)

    for i in 1:10
        QR(state, 1, 5, 9, 13)
        QR(state, 2, 6, 10, 14)
        QR(state, 3, 7, 11, 15)
        QR(state, 4, 8, 12, 16)

        QR(state, 1, 6, 11, 16)
        QR(state, 2, 7, 12, 13)
        QR(state, 3, 8, 9, 14)
        QR(state, 4, 5, 10, 15)
    end

    for i in 1:16
        state[i] += initialState[i]
    end
    return state
end